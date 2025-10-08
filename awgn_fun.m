function [y_copy,var] = awgn_fun(sig,reqSNR,ind,varargin)
%here what we want to do is only pickup the data part and the signal and
%add noise to that specific part, not to the LTF part of the signal

% creating the copy
y_copy = sig;

narginchk(2,6);

% Validate signal input
validateattributes(sig, {'numeric'}, ...
    {'nonempty'}, 'awgn', 'signal input');

dlarraySig = isa(sig,'dlarray');
gpuArraySig = isa(sig,'gpuArray');

% Formatted dlarrays are not supported
if dlarraySig && ~isempty(dims(sig))
    error(message('comm:checkinp:InvalidDlarrayFormat'));
end

numBatches = size(sig,3);

% Validate SNR input
validateattributes(reqSNR, {'numeric'}, ...
    {'real','vector','nonempty'}, 'awgn', 'SNR input');

vectorSNR = ~isscalar(reqSNR);

% Dlarray/gpuArray SIG do NOT support ndims(sig) > 3
if (dlarraySig || gpuArraySig || vectorSNR) & ndims(sig)>3
    error(message('comm:awgn:InvalidInputSize'));
end

if vectorSNR && length(reqSNR)~=numBatches
    error(message('comm:awgn:InvalidSNRSize',numBatches));
end

% Validate signal power
if nargin >= 4
    if strcmpi(varargin{1}, 'measured')
        if vectorSNR || dlarraySig || gpuArraySig
            error(message('comm:awgn:UnsupportedMeasuredSyntax'));
        end
        sigPower = sum(abs(sig(:)).^2)/numel(sig); % linear
    else
        validateattributes(varargin{1}, {'numeric'}, ...
            {'real','scalar','nonempty'}, 'awgn', 'signal power input');
        sigPower = varargin{1}; % linear or dB
    end
else
    sigPower = 1; % linear, default
end

% Validate state or power type
if nargin >= 5
    if comm.internal.utilities.isCharOrStringScalar(varargin{2}) && ...
            all(~strcmpi(varargin{2}, {'db','linear'}))
        error(message('comm:awgn:InvalidPowerType'));
    end

    isStream = ~isempty(varargin{2}) && ~comm.internal.utilities.isCharOrStringScalar(varargin{2});
    if isStream && (gpuArraySig||dlarraySig)
        error(message('comm:awgn:UnsupportedSeedStreamSyntax'));
    end

    if isStream && ~isa(varargin{2}, 'RandStream') % Random stream seed
        validateattributes(varargin{2}, {'double'}, ...
            {'real','scalar','nonnegative','integer','<',2^32}, ...
            'awgn', 'seed input');
    end
else % Use default stream & seed
    isStream = false;
end

% Validate power type
if nargin >= 5
    if comm.internal.utilities.isCharOrStringScalar(varargin{2}) % Type has been specified as the 4th input
        error(message('comm:awgn:InputAfterPowerType'));
    end
    if all(~strcmpi(varargin{3}, {'db','linear'}))
        error(message('comm:awgn:InvalidPowerType'));
    end
end

isLinearScale = ((nargin == 5) && ~isStream && strcmpi(varargin{3}, 'linear')) || ...
    ((nargin == 6) && strcmpi(varargin{4}, 'linear'));

% Cross-validation
if isLinearScale && (sigPower < 0)
    error(message('comm:awgn:InvalidSigPowerForLinearMode'));
end

if isLinearScale && any(reqSNR < 0)
    error(message('comm:awgn:InvalidSNRForLinearMode'));
end

if ~isLinearScale  % Convert signal power and SNR to linear scale
    if (nargin >= 4) && ~comm.internal.utilities.isCharOrStringScalar(varargin{1}) % User-specified signal power
        sigPower = 10^(sigPower/10);
    end
    reqSNR = 10.^(reqSNR./10);
end

if vectorSNR
    noisePower = reshape(sigPower./reqSNR,1,1,numBatches);
else
    noisePower = sigPower./reqSNR;
end

%now lets pick up only the data indices
data_ind = ind.EHTData ;
start = data_ind(1);
endind = data_ind(2);
if isStream
    if isa(varargin{3}, 'RandStream')
        stream = varargin{3};
    elseif isempty(coder.target)
        stream = RandStream('mt19937ar', 'Seed', varargin{3});
    else
        stream = coder.internal.RandStream('mt19937ar', 'Seed', varargin{3});
    end
    y_copy(start : endind) = sig(start : endind) + sqrt(noisePower).*randn(stream,size(sig(start : endind)),"like",sig);
else
    y_copy(start : endind) = sig(start : endind) + sqrt(noisePower).*randn(size(sig(start : endind)),"like",sig);
end

if nargout == 3
    var = reshape(noisePower, size(reqSNR));
end
