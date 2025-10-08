function [y , g] = ehtData_new(txPSDU,cfg,varargin)
%ehtData EHT Data field carrying the PSDU(s)
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   Y = ehtData(PSDU,CFG) generates the EHT Data field carrying the
%   PSDU(s) time-domain signal for the EHT transmission format.
%
%   Y is the time-domain EHT Data field signal. It is a complex matrix of
%   size Ns-by-Nt, where Ns represents the number of time-domain samples
%   and Nt represents the number of transmit antennas.
%
%   PSDU is the PHY service data unit input to the PHY. For single-user,
%   PSDU can be a double or int8 typed binary column vector of length
%   cfg.getPSDULength*8.
%
%   CFGEHT is the format configuration object of type <a href="matlab:help('wlanEHTMUConfig')">wlanEHTMUConfig</a> or 
%   <a href="matlab:help('wlanEHTTBConfig')">wlanEHTTBConfig</a>.
%
%   Y = ehtData(...,SCRAMINIT) optionally allows specification of the
%   scrambler initialization.
%
%   When not specified, SCRAMINIT defaults to a value of 93. When
%   specified, it is a double or int8-typed integer scalar or 1-by-Nu row
%   vector between 1 and 2047, inclusive, where Nu represents the number of
%   users.
%
%   Y = ehtData(...,SCRAMINIT,OSF) optionally allows specification of the
%   oversampling factor. When not specified, OSF defaults to 1.

%   Copyright 2022 The MathWorks, Inc.

%#codegen

cbw = wlan.internal.cbwStr2Num(cfg.ChannelBandwidth);
allocationInfo = ruInfo(cfg);
numRUs = allocationInfo.NumRUs;
ruSizes = allocationInfo.RUSizes;
ruIndices = allocationInfo.RUIndices;
numUsers = allocationInfo.NumUsers;
alpha = allocationInfo.PowerBoostFactorPerRU;
numSTSTotalRU = allocationInfo.NumSpaceTimeStreamsPerRU;
osf = 1;

% Default baseband rate
if nargin == 2
    scramInitBits = uint8(repmat([0; 0; 0; 0; 1; 0; 1; 1; 1; 0; 1],1,numUsers)); % Default is 93
else
    % Validate scrambler initialization input
    scramInitBits = varargin{1};
    if nargin == 4
       osf = varargin{2}; 
    end
end
[commonCodingParams,userCodingParams] = wlan.internal.ehtCodingParameters(cfg);

numTx = cfg.NumTransmitAntennas;
% Get OFDMA tone allocation constants for the RUs/MRUs

Nfft = 256*cbw/20;
ofdmGrid = complex(zeros(Nfft,commonCodingParams.NSYM,numTx));
isEHTMU = strcmp(packetFormat(cfg),'EHT-MU');
NUSIG = 2; % Number of OFDM symbols in U-SIG field

if isEHTMU
    sigInfo = wlan.internal.ehtSIGCodingInfo(cfg);
    NEHTSIG = sigInfo.NumSIGSymbols; % Number of OFDM symbols in EHT-SIG field
    z = 2+NUSIG+NEHTSIG; % Pilot symbol offset

    if ~iscell(txPSDU)
        txPSDUCell = {txPSDU};
    else
        txPSDUCell = txPSDU;
    end

    % Frequency and spatial mapping
    cardKr = coder.nullcopy(zeros(1,numRUs)); % Resource unit size for RUs/MRUs
    for j = 1:numRUs

        ruToneConstants = wlan.internal.heRUToneAllocationConstants(sum(ruSizes{j}));
        cardKr(j) = ruToneConstants.NST;
        % Get the index of the RU object containing the active RU properties
        iru = allocationInfo.RUNumbers(j);

        % Get current RU assignment details
        ru = cfg.RU{iru};
        numSTSTotal = allocationInfo.NumSpaceTimeStreamsPerRU(j);
        numUsers = allocationInfo.NumUsersPerRU(j);

        ruGrid = coder.nullcopy(complex(zeros(ruToneConstants.NST,commonCodingParams.NSYM,numSTSTotal)));

        % Get indices of data and pilots within RU (without Nulls) and active frequency indices.
        [ruMappingInd,k] = wlan.internal.ehtOccupiedSubcarrierIndices(cbw,ru.Size,ru.Index);

        % Encode, modulate, and map each user within the RU to space-time streams
        startSpaceTimeStream = 1;
        for i = 1:numUsers
            u = ru.UserNumbers(i);
            user = cfg.User{u};

            if user.APEPLength~=0
                if user.PostFECPaddingSource~=wlan.type.PostFECPaddingSource.userdefined
                    % Generate post-FEC padding bits
                    postFECPaddingBits = randomPostFECPadding(userCodingParams(u).NPADPostFEC,user.PostFECPaddingSource,user.PostFECPaddingSeed);
                else
                    % Repeat/cut/reshape user input bits
                    postFECPaddingBits = wlan.internal.parseInputBits(int8(user.PostFECPaddingBits),userCodingParams(u).NPADPostFEC);
                end

                % Generate STS for the user
                stsPerUser = wlan.internal.ehtSTSPerUser(txPSDUCell{u},scramInitBits(:,u),postFECPaddingBits,ru.Size,userCodingParams(u)); % Needs RU sizes and order

                % Determine the indices of space-time streams to map
                stsIdx = startSpaceTimeStream-1+(1:user.NumSpaceTimeStreams);
                % Map user data and pilots
                ruGrid(ruMappingInd.Data,:,stsIdx) = stsPerUser;
                n = 0:commonCodingParams.NSYM-1;
                ruGrid(ruMappingInd.Pilot,:,stsIdx) = wlan.internal.ehtPilots(ru.Size,user.NumSpaceTimeStreams,n,z); % Needs RU sizes and order
            end
            startSpaceTimeStream = startSpaceTimeStream+user.NumSpaceTimeStreams;
        end

        % Cyclic shift
        ruGrid = wlan.internal.heCyclicShift(ruGrid,cbw,k);

        % Spatial mapping
        txRUGrid = wlan.internal.spatialMap(ruGrid,ru.SpatialMapping,numTx,ru.SpatialMappingMatrix);

        % Calculate per RU scaling
        ruScalingFactor = alpha(j)/sqrt(numSTSTotalRU(j));

        % Map subcarriers to full FFT grid and scale
        ofdmGrid(k+Nfft/2+1,:,:) = txRUGrid*ruScalingFactor;
    end

    % Overall scaling factor
    allScalingFactor = Nfft/sqrt(sum(alpha.^2.*cardKr));
else
    ruToneConstants = wlan.internal.heRUToneAllocationConstants(sum(ruSizes{1}));
    cardKr = ruToneConstants.NST;

    NEHTSIG = 0;
    z = 2+NUSIG+NEHTSIG; % Pilot symbol offset

    if iscell(txPSDU)
        txPSDUMat = txPSDU{1};
    else
        txPSDUMat = txPSDU;
    end

    if cfg.PostFECPaddingSource~=wlan.type.PostFECPaddingSource.userdefined
        % Generate post-FEC padding bits
        postFECPaddingBits = randomPostFECPadding(userCodingParams.NPADPostFEC,cfg.PostFECPaddingSource,cfg.PostFECPaddingSeed);
    else
        % Repeat/cut/reshape user input bits
        postFECPaddingBits = reshape(...
        wlan.internal.parseInputBits(int8(cfg.PostFECPaddingBits),userCodingParams.NPADPostFEC), ...
        userCodingParams.NPADPostFEC,userCodingParams.mSTBC);
    end

    % Generate per-space-time stream symbols for user
    stsPerUser = wlan.internal.ehtSTSPerUser(txPSDUMat,scramInitBits,postFECPaddingBits,ruSizes{1},userCodingParams);

    % Map user data and pilots
    % Get indices of data and pilots within RU (without Nulls) and active frequency indices
    [ruMappingInd,k] = wlan.internal.ehtOccupiedSubcarrierIndices(cbw,ruSizes{1},ruIndices{1});
    ruGrid = coder.nullcopy(complex(zeros(cardKr,commonCodingParams.NSYM,cfg.NumSpaceTimeStreams)));
    ruGrid(ruMappingInd.Data,:,:) = stsPerUser(:,:,1:cfg.NumSpaceTimeStreams);
    n = 0:commonCodingParams.NSYM-1;
    ruGrid(ruMappingInd.Pilot,:,:) = wlan.internal.ehtPilots(ruSizes{1},cfg.NumSpaceTimeStreams,n,z);

    % Cyclic shift
    stsIdx = cfg.StartingSpaceTimeStream-1+(1:cfg.NumSpaceTimeStreams).'; 
    ruGrid = wlan.internal.heCyclicShift(ruGrid,cbw,k,stsIdx);

    % Spatial mapping
    Q = wlan.internal.ehtExtractRUFromSpatialMappingMatrix(cfg);
    txGrid = wlan.internal.spatialMap(ruGrid,cfg.SpatialMapping,numTx,Q);

    % Calculate per RU scaling
    ruScalingFactor = 1/sqrt(numSTSTotalRU);

    % Map subcarriers to full FFT grid and scale
    ofdmGrid(k+Nfft/2+1,:,:) = txGrid*ruScalingFactor;

    % Overall scaling factor
    allScalingFactor = Nfft/sqrt(cardKr);
end

% OFDM modulate and scale
switch cfg.GuardInterval
    case 0.8
        cpLen = 0.8*cbw;
    case 1.6
        cpLen = 1.6*cbw;
    otherwise % 3.2
        cpLen = 3.2*cbw;
end
[y , g] = wlanOFDMModulate_new(ofdmGrid,cpLen,osf);
y = y.*allScalingFactor;
g = g.*allScalingFactor;
end


function postFECPPadBits = randomPostFECPadding(NPADPostFEC,randomStream,seed)
    persistent stream

    if randomStream==wlan.type.PostFECPaddingSource.mt19937arwithseed
        if isempty(stream)
            if coder.target('MATLAB')
                stream = RandStream('mt19937ar');
            else
                stream = coder.internal.RandStream('mt19937ar');
            end
        end
        % Use same seed when generating each packet
        reset(stream,seed);
        postFECPPadBits = randi(stream,[0 1],NPADPostFEC,1,'int8');
    else
        % Use global stream
        postFECPPadBits = randi([0 1],NPADPostFEC,1,'int8');
    end
end