function [yout , postShift] = modulate_ofdm(gridIn,prmStr)
% COMM.INTERNAL.OFDM.MODULATE OFDM modulation
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   Y = comm.internal.ofdm.modulate(GRIDIN,PRMSTR) performs OFDM modulation
%   on the input GRIDIN, using the parameters specified in the structure
%   PRMSTR, and returns the result in Y. GRIDIN is the fully populated
%   3D input accounting for all data, null and pilot subcarriers.
%
%   PRMSTR must have the following fields:
%       FFTLength
%       CyclicPrefixLength
%       NumSymbols
%       NumTransmitAntennas
%   PRMSTR may have the following optional fields:
%       OversamplingFactor
%       NumBatchObs
%
%   See also comm.internal.ofdm.demodulate, ofdmdemod, ofdmmod.

%   Copyright 2017-2023 The MathWorks, Inc.

%#codegen

    fftLen = prmStr.FFTLength;
    cpLen  = prmStr.CyclicPrefixLength;
    numSym = prmStr.NumSymbols;
    numTx  = prmStr.NumTransmitAntennas;
    [osf, numBatchObs] = setup(prmStr);

    typeOut = cast(1i, 'like', gridIn);
    
    fftLenOSF = round(osf*fftLen);     % oversampled FFT length 
    cpLenOSF  = round(osf*cpLen);      % oversampled CP length(s)

    % Shift and IFFT
    zeroPaddedGrid = zeros(fftLenOSF,numSym,numTx,numBatchObs,'like',typeOut);
    if isreal(gridIn)
        zeroPaddedGrid(ceil((fftLenOSF-fftLen)/2)+(1:fftLen),:,:,:) = complex(gridIn,0);
    else
        zeroPaddedGrid(ceil((fftLenOSF-fftLen)/2)+(1:fftLen),:,:,:) = gridIn;
    end
    %% the postShift variable constains what we want, that contains the modulated symbols, in grid format.
    postShift = ifftshift(zeroPaddedGrid,1);
    postIFFT = osf*ifft(postShift,[],1); 
        
    % Append cyclic prefix
    if isscalar(cpLen) % same length
        % cpLen <= fftLen
        cpLen1OSF = cpLenOSF(1);   % oversampled CP length (first CP entry)
        postCP = postIFFT([end-cpLen1OSF+(1:cpLen1OSF),1:end],:,:,:);
        yout = reshape(postCP,[(fftLenOSF+cpLen1OSF)*numSym numTx numBatchObs]);

    else % different lengths per symbol
        
       yout = coder.nullcopy( ...
           zeros(fftLenOSF*numSym+sum(cpLenOSF),numTx,numBatchObs,'like',typeOut));
        for symIdx = 1:numSym
            cpLenOSFidx = cpLenOSF(symIdx);
            if cpLenOSFidx > fftLenOSF
                intSyms = floor(cpLenOSFidx/fftLenOSF);   % num full repeats of symbol
                newCPLenOSF = mod(cpLenOSFidx,fftLenOSF); % new cplen for symbol
                yout(fftLenOSF*(symIdx-1)+sum(cpLenOSF(1:symIdx-1)) + ...
                     (1:fftLenOSF+cpLenOSFidx),:,:) = ...
                    [reshape(postIFFT(end-newCPLenOSF+(1:newCPLenOSF), ...
                     symIdx,:,:), [newCPLenOSF numTx numBatchObs]); ...
                     reshape(repmat(postIFFT(:,symIdx,:,:),intSyms+1,1),[],numTx,numBatchObs)];
            else
                % Use reshape instead of squeeze in case of CP length ==1
                yout(fftLenOSF*(symIdx-1)+sum(cpLenOSF(1:symIdx-1)) + ...
                     (1:fftLenOSF+cpLenOSFidx),:,:) = ...
                    [reshape(postIFFT(end-cpLenOSFidx+(1:cpLenOSFidx), ...
                    symIdx,:,:), [cpLenOSFidx numTx numBatchObs]); ...
                    reshape(postIFFT(:,symIdx,:,:), [],numTx,numBatchObs)];
            end
        end

    end

end

function [osf, numBatchObs] = setup(prmStr)
    % check if OverSamplingFactor is specified, otherwise set it to 1
    if ~isfield(prmStr,'OversamplingFactor')
        osf = 1;
    else
        osf = prmStr.OversamplingFactor;
    end
    
    % Check if numBatchObs is specified, otherwise set it to 1
    if isfield(prmStr,'NumBatchObs')
        numBatchObs = prmStr.NumBatchObs;
    else
        numBatchObs = 1;
    end
end

% [EOF]
