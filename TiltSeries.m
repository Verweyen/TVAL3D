classdef TiltSeries < handle
    properties
        MaxMeasure 
        MinMeasure
        MeasurementFile
        AngleFile
        RawMeasurements
        Measurements
        TiltAngles
        NumAngles
        Dimensions
        TargetInteriorVolumes
    end
    methods
        function obj = TiltSeries(arg1, arg2)
        %%% TiltSeries::TiltSeries
        % Constructor for the TiltSeries object. Takes two arguments, one to 
        % specify the raw tilt series measurements, and the second to specify
        % the angles at which the tilt series was measured.
        %
        % TiltSeries(MeasurementFile,AngleFile)
        % TiltSeries(MeasurementTensor,AngleVector)
        obj.MeasurementFile = [];
        obj.AngleFile = [];
            if ischar(arg1)
                if ~exist(arg1,'file')
                    error('Specified tilt series measurement file <%s> does not exist.\n',arg1);
                end
                if ~exist(arg2,'file')
                    error('Specified tilt series measurement file <%s> does not exist.\n',arg2);
                end
                obj.MeasurementFile = arg1;
                obj.AngleFile = arg2;    
                obj.MaxMeasure = 255;
                obj.MinMeasure = 0;
                obj.TargetInteriorVolumes = [];
                obj.ReloadData;
            else
                if nargin < 2
                    error('Specify constructor as TiltSeries(Measurements,Angles)\n.');
                end
                if isscalar(arg1)
                    % In this case, we have been given a side-length for
                    % each of the observations and we should initialize all
                    % measurements to the zero case.
                    obj.RawMeasurements = zeros(arg1,arg1,length(arg2));
                else
                    obj.RawMeasurements = arg1;
                end
                obj.MaxMeasure = max(arg1(:));                
                obj.MinMeasure = max(arg1(:));                
                obj.TargetInteriorVolumes = [];
                obj.TiltAngles = arg2;
                obj.Dimensions = size(obj.RawMeasurements);
            end
        end

        function ReloadData(obj)
        %%% TiltSeries::ReloadData
        % Given that the measurement and angle files were specified in 
        % the constructor, this function will reload its internal storage
        % with the data on disk.
            if isempty(obj.MeasurementFile) || isempty(obj.AngleFile)
                error('Cannot reload tilt series from disk if file locations not specified.\n');
            else
                obj.TiltAngles = load(obj.AngleFile);
                obj.NumAngles = length(obj.TiltAngles);   
                % Loop through all of the tilt observations and save the raw image
                for angidx = 1:obj.NumAngles
                    obj.RawMeasurements(:,:,angidx) = double(imread(obj.MeasurementFile,'Index',angidx));
                end
                obj.Dimensions = size(obj.RawMeasurements);
                obj.Measurements = obj.RawMeasurements;
            end
            obj.MaxMeasure = 255;
            obj.MinMeasure = 0;
            obj.TargetInteriorVolumes = 1:1:obj.Dimensions(1);
        end

        function Reset(obj)
        %%% TiltSeries::Reset
        % Undo any modificaitons to the Measurements by copying RawMeasurements
        % back into the Measurements field
            obj.Measurements = obj.RawMeasurements;
            obj.MaxMeasure = max(obj.Measurements(:));
            obj.MinMeasure = min(obj.Measurements(:));
            obj.TargetInteriorVolumes = [];
            obj.Dimensions = size(obj.Measurements);
        end

        function Resize(obj,m,n,ResizeFilter)
        %%% TiltSeries::Resize
        % Resize the measurements according to the dimensions specified in
        % row-column order.
            tmpMeasurements = zeros(m,n,obj.NumAngles);
            for angidx=1:obj.NumAngles
                angleimage = squeeze(obj.Measurements(:,:,angidx));
                if nargin > 3
                    tmpMeasurements(:,:,angidx) = imresize(angleimage,[m,n],ResizeFilter);
                else
                    tmpMeasurements(:,:,angidx) = imresize(angleimage,[m,n]);
                end
            end
            obj.Measurements = tmpMeasurements;
            obj.Dimensions = [m,n,obj.NumAngles];
            obj.TargetInteriorVolumes = 1:1:obj.Dimensions(1);
        end

        function Normalize(obj)
        %%% TiltSeries::Normalize
        % Normalize the tilt series measurements to be in the range [0,1]
            maxMeasure = max(obj.Measurements(:));
            minMeasure = min(obj.Measurements(:));
            obj.Measurements = (obj.Measurements - minMeasure)./(maxMeasure - minMeasure);
            obj.MaxMeasure = 1;
            obj.MinMeasure = 0;
        end

        function X = FaceSlice(obj,SliceType,index)
        %%% TiltSeries::FaceSlice
        % Create an en-face slice through the tilt series measurements according to
        % the specified index
            SliceType = lower(SliceType);
            switch SliceType
                case {'h','horz','horizontal'}
                    X = squeeze(obj.Measurements(index,:,:));
                case {'v','vert','vertical'}
                    X = squeeze(obj.Measurements(:,index,:));
                case {'a','ang','angle'}
                    X = squeeze(obj.Measurements(:,:,index));
            end
        end

        function X = AtAngle(obj,fAngle)
        %%% TiltSeries::AtAngle
        % Return the measurement corresponding to the specified angle
            AngleDiffs = abs(obj.TiltAngles - fAngle);
            FoundLocations = AngleDiffs < 1e-8;
            if sum(FoundLocations)
                % If we find a negligibly close angle
                AngleIdx = find(FoundLocations);
                X = squeeze(obj.Measurements(:,:,AngleIdx));
            else
                error('No measurement in tilt series at angle <%0.4f deg>\n',fAngle)
            end
        end

        function X = Sinogram(obj,Line)
        %%% TiltSeries::Sinogram
        % Return the sinogram at a particular line. This is really just a 
        % wrapper around the FaceSlice function
            X = obj.FaceSlice('h',Line);
        end

        function X = PaddedSinogram(obj,Line)
        %%% TiltSeries::Padded Sinogram
        % Return a sinogram with additional padding on the line-measurements
        % in order to cooperate with Matlab's radon & iradon.            
            Sino = obj.Sinogram(Line);
            LineMeasurements = size(Sino,1);
            PaddedLineMeasurements = size(radon(zeros(LineMeasurements),0),1);

            top_pad = floor((PaddedLineMeasurements-LineMeasurements)./2);
            bot_pad = top_pad;
            if (top_pad + bot_pad + LineMeasurements) ~= PaddedLineMeasurements
                bot_pad = bot_pad + 1;
            end
            c = size(Sino,2);
            X = [ zeros(top_pad,c); ...
                  Sino; ...
                  zeros(bot_pad,c)];
        end

        function X = ctranspose(obj)
        %%% TiltSeries::ctranspose
        % Redefine the transpose operator for the object so that we
        % can directly get the back-projected volume from the the
        % tilt series.
        %
        % For right now, we will be using the back-projection provided
        % by Matlab's iradon, operating in the un-filtered mode.
            X = zeros(obj.Dimensions(1),obj.Dimensions(2));
            for i=1:length(obj.TargetInteriorVolumes)
                VolumeIdx = obj.TargetInteriorVolumes(i);
                X(:,:,i) = iradon(obj.FaceSlice('h',VolumeIdx),obj.TiltAngles,'linear',...
                                                                              'None',...
                                                                               1,...
                                                                               obj.Dimensions(1));
            end
        end

    end
end