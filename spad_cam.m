classdef spad_cam < handle
    %SPAD_CAM class to control SPAD512 camera from PiImaging
    %   Documentation of camera and commands here: https://piimaging.com/doc-cspad
    %

    properties
        host;       % Should be 'localhost' or '127.0.0.1'
        port;       % Defined in settings in CSPAD software
        connection; % Object used to send commands to SPAD interface server

        cleanupObj = onCleanup(@() obj.clearTcpClient()); % Clean up and close connection in case camera object is deleted
        default_timeout_seconds = 10;  % How long to wait for server response

        save_path;  % Where data will be saved. In save_path/data/intensity_images/acqXXXXX
        msg;        % Byte-wise message response from server
        msgstr;     % Decoded text of msg
        verbose = false; % If you want to print all return messages or not

        frame_cap = 100000; % Set a maximum frame cap to send for a single acquisition

        fan_speed = 1;                  % Fan speed, 0 = off, 1 = on.
        bits = 1;                       % Bit depth of images, either 1,4,6,7,8,9,10,11,12
        measurement_time = 1;           % Frame exposure time. For 1 bit, it is in micro-seconds, for other, in milliseconds
        frames = 10000;                 % Number of frames to capture
        frame_trigger = 0;              % If you are using an external trigger
        overlap = 0;                    % If you want to overlap exposure and readout
        sparse = 0;                     % If you want to save data in sparse format. Only available for 1-bit mode
        stream = 0;                     % If you want to stream the data directly to the tcp connection (1) or save to disk (0)
        image_width = '';               % *OPTIONAL* Define the image width, from 1 to 512. If not defined, the camera captures the entire frame (512x512)
        image_width_prev = '';          % Stores image width defined in previous capture
        additional_trigger_delay = '';  % *OPTIONAL* Define a delay between frame trigger and exposure opening

        hot_pixels = [];                % Array defining hot pixels. 0 = ok pixels, 1 = hot pixel
        hot_pixels_lin = [];            % Linearly-reshaped hot-pixel array


    end

    properties (Dependent, SetAccess = private)
        im_size
        im_size_prev
    end


    methods(Static)

        function keepVec = make_keep_mask_fast(hot_pixels_lin, im_size, pixelOrder)
            %MAKE_KEEP_MASK_FAST Convert hot-pixel input into a keep mask.
            %
            % hot_pixels_lin supports:
            %   []                         -> keep all pixels
            %   image-sized mask           -> 1 keep, 0 discard
            %   vector of length prod(...) -> keep mask in the chosen pixel order
            %   shorter vector             -> MATLAB linear indices to discard
            %
            % pixelOrder:
            %   'row-major' -> stream order matching original packed-bit reader
            %   'matlab'    -> MATLAB column-major image linear order

            im_size = im_size(:).';
            nPix = prod(im_size);

            if nargin < 3 || isempty(pixelOrder)
                pixelOrder = 'row-major';
            end

            pixelOrder = normalize_pixel_order_fast(pixelOrder);

            if isempty(hot_pixels_lin)
                keepVec = true(nPix, 1);
                return;
            end

            if isequal(size(hot_pixels_lin), im_size)
                keepImage = hot_pixels_lin ~= 0;

                switch pixelOrder
                    case 'row-major'
                        keepVec = reshape(keepImage.', [], 1);
                    case 'matlab'
                        keepVec = keepImage(:);
                end

                return;
            end

            if numel(hot_pixels_lin) == nPix
                keepVec = hot_pixels_lin(:) ~= 0;
                return;
            end

            badIdx = double(hot_pixels_lin(:));

            bad = badIdx < 1 | badIdx > nPix | badIdx ~= floor(badIdx);
            if any(bad)
                error('hot_pixels_lin contains invalid MATLAB linear indices.');
            end

            keepImage = true(im_size);
            keepImage(badIdx) = false;

            switch pixelOrder
                case 'row-major'
                    keepVec = reshape(keepImage.', [], 1);
                case 'matlab'
                    keepVec = keepImage(:);
            end
        end

        function pixelOrder = normalize_pixel_order_fast(pixelOrder)
            %NORMALIZE_PIXEL_ORDER_FAST Normalize pixel-order option.

            key = lower(strtrim(char(pixelOrder)));
            key = strrep(key, '_', '-');
            key = strrep(key, ' ', '-');

            switch key
                case {'row-major', 'rowmajor', 'stream', 'c', 'corder', 'c-order'}
                    pixelOrder = 'row-major';

                case {'matlab', 'column-major', 'columnmajor', 'col-major', 'column', 'fortran'}
                    pixelOrder = 'matlab';

                otherwise
                    error('Unknown PixelOrder: %s', char(pixelOrder));
            end
        end
    end
    methods

        %% Basic camera operation
        function obj = spad_cam(host,port)
            %SPAD_CAM Construct an instance of this class
            %   Detailed explanation goes here
            obj.host = host;
            obj.port = port;
            obj.connection = tcpclient(obj.host, obj.port, ...
                'Timeout', 30, ...
                'ConnectTimeout', 10, ...
                'EnableTransferDelay', false);
            obj.readResponse;
            disp(obj.msgstr);

            obj.setFanSpeed;

        end

        function clearTcpClient(obj)
            % Delete the TCP client when the script finishes or encounters an error.
            if ~isempty(obj.connection) && isvalid(obj.connection)
                delete(obj.connection);
            end
        end

        function sendCommand(obj, command)
            %
            %   Function to send commands to SPAD camera.
            %
            %   command is a string composed of an initial character/word
            %   defining the specific command, followed by any parameters
            %   to go with that command.
            %
            %   See docs for more details on commands and command
            %   structure: https://piimaging.com/doc-cspad


            % Append an actual newline byte.
            %
            % Using 'S,2\n' directly in a MATLAB character vector would normally
            % send the two literal characters "\" and "n". sprintf converts \n
            % into an actual newline character.
            command = uint8(char(command));
            data = uint8(sprintf('%s\n', command));

            write(obj.connection, data, 'uint8');
        end

        function readResponse(obj)
            % Wait for the first response byte.

            timeoutSeconds = obj.default_timeout_seconds;

            startTime = tic;
            while obj.connection.NumBytesAvailable == 0
                if toc(startTime) > timeoutSeconds
                    error('TCP:ResponseTimeout', ...
                        'No server response was received within %.1f seconds.', ...
                        timeoutSeconds);
                end

                pause(0.05);
            end

            % Wait briefly for the rest of the response to arrive. This avoids
            % reading only the first TCP packet of a multi-packet message.
            previousCount = -1;
            stableChecks = 0;

            while stableChecks < 3
                currentCount = obj.connection.NumBytesAvailable;

                if currentCount == previousCount
                    stableChecks = stableChecks + 1;
                else
                    stableChecks = 0;
                    previousCount = currentCount;
                end

                pause(0.05);
            end

            % Read all currently available bytes.
            obj.msg = read(obj.connection, obj.connection.NumBytesAvailable, 'uint8');

            % Convert bytes to MATLAB text.
            obj.msgstr = native2unicode(obj.msg, 'UTF-8');

            % Remove trailing CR/LF characters from the displayed response.
            obj.msgstr = strtrim(obj.msgstr);
        end

        function setFanSpeed(obj,varargin)
            if length(varargin)>0
                obj.fan_speed = varargin{1};
            end
            fan_string = ['S,' num2str(obj.fan_speed)];
            obj.sendCommand(fan_string);
            obj.readResponse;
            disp(obj.msgstr);
        end

        function setSaveFolder(obj,save_path)
            obj.save_path = save_path;
            change_path_string = ['D,' obj.save_path];

            obj.sendCommand(change_path_string);
            obj.readResponse;

            disp('Set saving folder to: ')
            disp(obj.msgstr);
        end

        function imageAcquisitionSingle(obj,n_frames)
            measurement_string = ['I,' ...
                num2str(obj.bits) ',' ...
                num2str(obj.measurement_time) ',' ...
                num2str(n_frames) ',' ...
                num2str(obj.frame_trigger) ',' ...
                num2str(obj.overlap) ',' ...
                num2str(obj.sparse) ',' ...
                num2str(obj.stream) ',' ...
                num2str(obj.image_width) ',' ...
                num2str(obj.additional_trigger_delay)
                ];

            disp('Sending acquisition...')
            obj.sendCommand(uint8(char(measurement_string)))
            obj.readResponse();

            disp('Finished acquisition.')

            if obj.verbose
                disp('Return message: ')
                disp(obj.msgstr);
            end

            obj.image_width_prev = obj.image_width;

        end

        function imageAcquisition(obj)
            % Take care of when the number of desired frames is greater
            % than the maximum that the camera can capture with one command
            % (130,000)

            if obj.frames<=obj.frame_cap
                obj.imageAcquisitionSingle(obj.frames)
            else

                n_sets = ceil(obj.frames/obj.frame_cap);

                disp(['Number of frames too large. Splitting into multiple acquisitions. Number of acquisitions = ' num2str(n_sets)])
                frames_captured = 0;

                for set = 1:n_sets
                    disp([' ______ Acqusition ' num2str(set) ' of ' num2str(n_sets) ' ______ '])
                    if set<n_sets
                        n_frames = obj.frame_cap;
                    else
                        n_frames = obj.frames-frames_captured;
                    end
                    obj.imageAcquisitionSingle(n_frames)
                    frames_captured = frames_captured+n_frames;
                    pause(2)
                end
            end

            obj.image_width_prev = obj.image_width;
        end
        
        %% Basic data reading 

        %% Correlations measurements 

        %% 
        
        %% Auto-updating of useful internal parameters
        function value = get.im_size(obj) 
            value = [obj.image_width 512];
        end

        function value = get.im_size_prev(obj)
            value = [obj.image_width_prev 512];
        end
   
    
    end
end