% =========================================================================   
% (c) 2018 Ronald Nissel, ronald.nissel@gmail.com
% =========================================================================      
% This script allows to reproduce Figure 7 of "FBMC-OQAM in Doubly-
% Selective Channels: A New Perspective on MMSE Equalization", R.Nissel,
% M.Rupp, R.Marsalek, IEEE SPAWC 2017. In particular it simulates the 2x2
% MIMO BER of OFDM and FBMC in doubly-selective channels and compares 
% different MMSE equalization methods as well as a low-complexity 
% interference cancellation scheme.

clear; close all;

%% Simulation Parameters
M_SNR_dB            = [10:2:40];           % SNR in dB
NrRepetitions       = 10;                  % Number of Monte Carlo repetitions (channel realizations) over which we average                  
NrRepetitions2      = 5000;                % Number of Monte Carlo repetitions for the doubly-flat reference curve

% Channel
Velocity_kmh        = 500;                 % Velocity in km/h. Note that [mph]*1.6=[kmh] and [m/s]*3.6=[kmh]
PowerDelayProfile   = 'VehicularA';        % Channel model, either string or vector: 'Flat', 'AWGN', 'PedestrianA', 'PedestrianB', 'VehicularA', 'VehicularB', 'ExtendedPedestrianA', 'ExtendedPedestrianB', or 'TDL-A_xxns','TDL-B_xxns','TDL-C_xxns' (with xx the RMS delay spread in ns, e.g. 'TDL-A_30ns'), or [1 0 0.2] (Self-defined power delay profile which depends on the sampling rate) 

% FBMC and OFDM
L                   = 24;                  % Number of subcarriers, 24 to keep the simulation time short. Moreover, a higher number does not give significant other results. 
F                   = 15e3;                % Subcarrier spacing in Hz, 15kHz, same as in LTE
SamplingRate        = F*14*14;             % Sampling rate in Hz. Must be a multiple of the subcarrier spacing. 14 because of the CP in OFDM. F*14*14 because the sampling rate should matche approximately the predefined channel delay taps (Vehicular A)
QAM_ModulationOrder = 256;                 % Modulation order, 4, 16, 64, 256, 1024, ...

NrIterationsInterferenceCancellation = 25; % Number of iterations for the simple interference cancellation scheme

%% Parameters in the paper
% NrRepetitions = 2000;
% NrRepetitions2 = 50000;

%% Additional Stuff
addpath('./Theory');
addpath('./EqualizerFunctions');
warning('off','MATLAB:nearlySingularMatrix'); % matrix inversion is near singular but still works. The error could be mitigated by using pinv(), which is based on a SVD. However, this takes longer so that I do not use it here.


%% FBMC Object
FBMC = Modulation.FBMC(...
    L,...                               % Number subcarriers
    30,...                              % Number FBMC symbols
    F,...                               % Subcarrier spacing (Hz)
    SamplingRate,...                    % Sampling rate (Samples/s)
    0,...                               % Intermediate frequency first subcarrier (Hz)
    false,...                           % Transmit real valued signal
    'Hermite-OQAM',...                  % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM, 
    8, ...                              % Overlapping factor (also determines oversampling in the frequency domain)
    0, ...                              % Initial phase shift
    true ...                            % Polyphase implementation
    );

%% OFDM Object (Add zeroes to the OFDM signal so that it fits the FBMC signal)
ZeroGuardTimeLength = ((FBMC.Nr.SamplesTotal-(round((1/15e3/14)*SamplingRate)+round(SamplingRate/15e3))*14)/2)/SamplingRate;
OFDM = Modulation.OFDM(...
    L,...                           % Number Subcarriers
    14,...                          % Number OFDM Symbols
    F,...                           % Subcarrier spacing (Hz)
    SamplingRate,...                % Sampling rate (Samples/s)
    0,...                           % Intermediate frequency first subcarrier (Hz)
    false,...                       % Transmitreal valued signal
    1/15e3/14, ...                  % Cyclic prefix length (s)
    ZeroGuardTimeLength ...         % Zero guard length (s)
    );

%% Check Number of Samples
if  OFDM.Nr.SamplesTotal~=FBMC.Nr.SamplesTotal
   error('Total number of samples must be the same for OFDM and FBMC.');
end
N = OFDM.Nr.SamplesTotal;


%% Channel Model Object
ChannelModel = Channel.FastFading(...
    SamplingRate,...                     % Sampling rate (Samples/s)
    PowerDelayProfile,...                % Power delay profile, either string or vector: 'Flat', 'AWGN', 'PedestrianA', 'PedestrianB', 'VehicularA', 'VehicularB', 'ExtendedPedestrianA', 'ExtendedPedestrianB', or 'TDL-A_xxns','TDL-B_xxns','TDL-C_xxns' (with xx the RMS delay spread in ns, e.g. 'TDL-A_30ns'), or [1 0 0.2] (Self-defined power delay profile which depends on the sampling rate) 
    N,...                                % Number of total samples
    Velocity_kmh/3.6*2.5e9/2.998e8,...   % Maximum Doppler shift: Velocity_kmh/3.6*CarrierFrequency/2.998e8  
    'Jakes',...                          % Which Doppler model: 'Jakes', 'Uniform', 'Discrete-Jakes', 'Discrete-Uniform'. For "Discrete-", we assume a discrete Doppler spectrum to improve the simulation time. This only works accuratly if the number of samples and the velocity is sufficiently large                                       
    200, ...                             % Number of paths for the WSSUS process. Only relevant for a 'Jakes' and 'Uniform' Doppler spectrum                                                 
    2,...                                % Number of transmit antennas
    2,...                                % Number of receive antennas
    true ...                             % Gives a warning if the predefined delay taps of the channel do not fit the sampling rate. This is usually not much of a problem if they are approximatly the same.
    );

%% PAM and QAM Object
PAM = Modulation.SignalConstellation(sqrt(QAM_ModulationOrder),'PAM');
QAM = Modulation.SignalConstellation(QAM_ModulationOrder,'QAM');


%% Precalculate Transmit and Receive Matrices
G_FBMC = FBMC.GetTXMatrix;
Q_FBMC = (FBMC.GetRXMatrix)'; % In FBMC the TX and RX matrices are the same (except a scaling factor)

G_OFDM = OFDM.GetTXMatrix;
Q_OFDM = (OFDM.GetRXMatrix)';

CorrelationMatrixNoise_FBMC = Q_FBMC'*Q_FBMC;
CorrelationMatrixNoise_OFDM = Q_OFDM'*Q_OFDM;

CorrelationMatrixNoise_OFDM_MIMO = [CorrelationMatrixNoise_OFDM,zeros(size(CorrelationMatrixNoise_OFDM));...
                                    zeros(size(CorrelationMatrixNoise_OFDM)),CorrelationMatrixNoise_OFDM];

CorrTemp = [Q_FBMC'*Q_FBMC,zeros(size(Q_FBMC,2));...
            zeros(size(Q_FBMC,2)),Q_FBMC'*Q_FBMC];
CorrelationMatrixNoise_FBMC_RealImag_MIMO =  1/2*[real(CorrTemp),-imag(CorrTemp);imag(CorrTemp),real(CorrTemp)];       
                                        
%% Preallocate for Parfor
BER_1TapMMSE_FBMC = nan(length(M_SNR_dB),NrRepetitions);
BER_1TapMMSE_OFDM = nan(length(M_SNR_dB),NrRepetitions);

BER_FullMMSE_FBMC = nan(length(M_SNR_dB),NrRepetitions);
BER_FullMMSE_OFDM = nan(length(M_SNR_dB),NrRepetitions);

BER_9TapTimeMMSE_FBMC = nan(length(M_SNR_dB),NrRepetitions);
BER_9TapMMSE_FBMC = nan(length(M_SNR_dB),NrRepetitions);

BER_InterferenceCancellation_OFDM = nan(length(M_SNR_dB),NrRepetitions,NrIterationsInterferenceCancellation);
BER_InterferenceCancellation_FBMC = nan(length(M_SNR_dB),NrRepetitions,NrIterationsInterferenceCancellation);

tic    
for i_rep = 1:NrRepetitions
    %% Update Channel
    ChannelModel.NewRealization;
    
    %% Binary Data
    BinaryDataStream_FBMC_Antenna1 = randi([0 1],FBMC.Nr.Subcarriers*FBMC.Nr.MCSymbols*log2(PAM.ModulationOrder),1);
    BinaryDataStream_FBMC_Antenna2 = randi([0 1],FBMC.Nr.Subcarriers*FBMC.Nr.MCSymbols*log2(PAM.ModulationOrder),1);

    BinaryDataStream_OFDM_Antenna1 = randi([0 1],OFDM.Nr.Subcarriers*OFDM.Nr.MCSymbols*log2(QAM.ModulationOrder),1);
    BinaryDataStream_OFDM_Antenna2 = randi([0 1],OFDM.Nr.Subcarriers*OFDM.Nr.MCSymbols*log2(QAM.ModulationOrder),1);
 
    %% Transmitted Data Symbols (Map bin to symbol)
    x_FBMC_Antenna1 = reshape(PAM.Bit2Symbol(BinaryDataStream_FBMC_Antenna1),FBMC.Nr.Subcarriers,FBMC.Nr.MCSymbols);
    x_FBMC_Antenna2 = reshape(PAM.Bit2Symbol(BinaryDataStream_FBMC_Antenna2),FBMC.Nr.Subcarriers,FBMC.Nr.MCSymbols);
    
    x_OFDM_Antenna1 = reshape(QAM.Bit2Symbol(BinaryDataStream_OFDM_Antenna1),OFDM.Nr.Subcarriers,OFDM.Nr.MCSymbols);
    x_OFDM_Antenna2 = reshape(QAM.Bit2Symbol(BinaryDataStream_OFDM_Antenna2),OFDM.Nr.Subcarriers,OFDM.Nr.MCSymbols);
        
    %% Transmitted FBMC Signal (time domain).
    s_FBMC_Antenna1 = G_FBMC*x_FBMC_Antenna1(:)/sqrt(2); % Same as "FBMC.Modulation(x_FBMX)" which is computationally more efficient. But G_FBMC is consistent to our paper.
    s_FBMC_Antenna2 = G_FBMC*x_FBMC_Antenna2(:)/sqrt(2);
    
    s_OFDM_Antenna1 = G_OFDM*x_OFDM_Antenna1(:)/sqrt(2);
    s_OFDM_Antenna2 = G_OFDM*x_OFDM_Antenna2(:)/sqrt(2);
    

    %% Channel
    ConvolutionMatrix11 = ChannelModel.GetConvolutionMatrix{1,1};
    ConvolutionMatrix12 = ChannelModel.GetConvolutionMatrix{1,2};
    ConvolutionMatrix21 = ChannelModel.GetConvolutionMatrix{2,1};
    ConvolutionMatrix22 = ChannelModel.GetConvolutionMatrix{2,2};
  
    r_FBMC_noNoise_Antenna1 = ConvolutionMatrix11*s_FBMC_Antenna1 + ConvolutionMatrix12*s_FBMC_Antenna2;
    r_FBMC_noNoise_Antenna2 = ConvolutionMatrix21*s_FBMC_Antenna1 + ConvolutionMatrix22*s_FBMC_Antenna2;
    
    r_OFDM_noNoise_Antenna1 = ConvolutionMatrix11*s_OFDM_Antenna1 + ConvolutionMatrix12*s_OFDM_Antenna2;
    r_OFDM_noNoise_Antenna2 = ConvolutionMatrix21*s_OFDM_Antenna1 + ConvolutionMatrix22*s_OFDM_Antenna2;
    
    %% For MMSE Equalization
    D_FBMC_MIMO = [Q_FBMC'*ConvolutionMatrix11*G_FBMC,Q_FBMC'*ConvolutionMatrix12*G_FBMC;...
                   Q_FBMC'*ConvolutionMatrix21*G_FBMC,Q_FBMC'*ConvolutionMatrix22*G_FBMC]/sqrt(2);
    D_OFDM_MIMO = [Q_OFDM'*ConvolutionMatrix11*G_OFDM,Q_OFDM'*ConvolutionMatrix12*G_OFDM;...
                   Q_OFDM'*ConvolutionMatrix21*G_OFDM,Q_OFDM'*ConvolutionMatrix22*G_OFDM]/sqrt(2);  

    %% Channel for one tap equalizer
    h_FBMC_MIMO = diag(D_FBMC_MIMO);
    h_OFDM_MIMO = diag(D_OFDM_MIMO);
    
    % parfor i_SNR = 1:length(M_SNR_dB)      
    for i_SNR = 1:length(M_SNR_dB)
        %% Add Noise
        SNR_dB = M_SNR_dB(i_SNR);
        Pn = 10^(-SNR_dB/10);
        Pn_time = SamplingRate/(F*L)*10^(-SNR_dB/10);
        noise_Antenna1 = sqrt(Pn_time/2)*(randn(N,1)+1j*randn(N,1));
        noise_Antenna2 = sqrt(Pn_time/2)*(randn(N,1)+1j*randn(N,1));

        r_FBMC_Antenna1 = r_FBMC_noNoise_Antenna1 + noise_Antenna1; 
        r_FBMC_Antenna2 = r_FBMC_noNoise_Antenna2 + noise_Antenna2; 

        r_OFDM_Antenna1 = r_OFDM_noNoise_Antenna1 + noise_Antenna1;
        r_OFDM_Antenna2 = r_OFDM_noNoise_Antenna2 + noise_Antenna2;

        %% Demodulate FBMC signal
        y_FBMC_Antenna1 = Q_FBMC'*r_FBMC_Antenna1; % Same as "FBMC.Demodulation(r_FBMC)" 
        y_FBMC_Antenna2 = Q_FBMC'*r_FBMC_Antenna2; % Same as "FBMC.Demodulation(r_FBMC)" 
        y_FBMC = [y_FBMC_Antenna1;y_FBMC_Antenna2];

        y_OFDM_Antenna1 = Q_OFDM'*r_OFDM_Antenna1;
        y_OFDM_Antenna2 = Q_OFDM'*r_OFDM_Antenna2;
        y_OFDM = [y_OFDM_Antenna1;y_OFDM_Antenna2];

        %%  MMSE One Tap Per Antenna Equalizer = 4Taps (interference due to doubly-selectivity is not included here)
        x_est_1TapMMSE_FBMC = nan(FBMC.Nr.Subcarriers*FBMC.Nr.MCSymbols,1);
        for i_lk = 1:FBMC.Nr.Subcarriers*FBMC.Nr.MCSymbols
            IndexFBMC = zeros(1,FBMC.Nr.Subcarriers*FBMC.Nr.MCSymbols*2);
            IndexFBMC([i_lk i_lk+FBMC.Nr.Subcarriers*FBMC.Nr.MCSymbols])=1;
            H_FBMC = D_FBMC_MIMO(IndexFBMC==1,IndexFBMC==1);
            x_est_1TapMMSE_FBMC(IndexFBMC==1) = real(H_FBMC'*(H_FBMC*H_FBMC' + eye(2)*Pn)^-1*y_FBMC(IndexFBMC==1));     
        end

        x_est_1TapMMSE_OFDM = nan(OFDM.Nr.Subcarriers*OFDM.Nr.MCSymbols,1);
        for i_lk = 1:OFDM.Nr.Subcarriers*OFDM.Nr.MCSymbols
            IndexOFDM = zeros(1,OFDM.Nr.Subcarriers*OFDM.Nr.MCSymbols*2);
            IndexOFDM([i_lk i_lk+OFDM.Nr.Subcarriers*OFDM.Nr.MCSymbols])=1;
            H_OFDM = D_OFDM_MIMO(IndexOFDM==1,IndexOFDM==1);     
            x_est_1TapMMSE_OFDM(IndexOFDM==1) = H_OFDM'*(H_OFDM*H_OFDM'+eye(2)*Pn)^-1*y_OFDM(IndexOFDM==1);     
        end

        %% Low-Complexity Interference Cancellation (without sorting etc.)
        x_est_FBMC_Temp = x_est_1TapMMSE_FBMC; % initialize with one tap estimates    
        for i_iterationFBMC = 1:NrIterationsInterferenceCancellation     
            x_est_FBMC_Temp2 = nan(size(x_est_FBMC_Temp));
            for i_lk = 1:FBMC.Nr.Subcarriers*FBMC.Nr.MCSymbols
                IndexFBMC = zeros(1,FBMC.Nr.Subcarriers*FBMC.Nr.MCSymbols*2);
                IndexFBMC([i_lk i_lk+FBMC.Nr.Subcarriers*FBMC.Nr.MCSymbols])=1;
                H_FBMC = D_FBMC_MIMO(IndexFBMC==1,IndexFBMC==1);
                D_Interference = D_FBMC_MIMO(IndexFBMC==1,IndexFBMC==0);
                y_FBMC_InterferenceCancellation_Temp = y_FBMC(IndexFBMC==1)-D_Interference*x_est_FBMC_Temp(IndexFBMC==0);
                x_est = real(H_FBMC'*(H_FBMC*H_FBMC' + eye(2)*Pn)^-1*y_FBMC_InterferenceCancellation_Temp);
                x_est_FBMC_Temp2(IndexFBMC==1) =  x_est;
            end
            x_est_FBMC_Temp = PAM.SymbolQuantization(x_est_FBMC_Temp2);        
            DetectedBitStream_FBMC_Temp = PAM.Symbol2Bit(x_est_FBMC_Temp);        
            BER_InterferenceCancellation_FBMC(i_SNR,i_rep,i_iterationFBMC) = mean([BinaryDataStream_FBMC_Antenna1;BinaryDataStream_FBMC_Antenna2]~=DetectedBitStream_FBMC_Temp);
        end

        x_est_OFDM_Temp = x_est_1TapMMSE_OFDM; % initialize with one tap estimates    
        for i_iterationOFDM = 1:NrIterationsInterferenceCancellation     
            x_est_OFDM_Temp2 = nan(size(x_est_OFDM_Temp));
            for i_lk = 1:OFDM.Nr.Subcarriers*OFDM.Nr.MCSymbols
                IndexOFDM = zeros(1,OFDM.Nr.Subcarriers*OFDM.Nr.MCSymbols*2);
                IndexOFDM([i_lk i_lk+OFDM.Nr.Subcarriers*OFDM.Nr.MCSymbols])=1;
                H_OFDM = D_OFDM_MIMO(IndexOFDM==1,IndexOFDM==1);
                D_Interference = D_OFDM_MIMO(IndexOFDM==1,IndexOFDM==0);
                y_OFDM_InterferenceCancellation_Temp = y_OFDM(IndexOFDM==1)-D_Interference*x_est_OFDM_Temp(IndexOFDM==0);
                x_est = H_OFDM'*(H_OFDM*H_OFDM' + eye(2)*Pn)^-1*y_OFDM_InterferenceCancellation_Temp;
                x_est_OFDM_Temp2(IndexOFDM==1) =  x_est;
            end
            x_est_OFDM_Temp = QAM.SymbolQuantization(x_est_OFDM_Temp2);        
            DetectedBitStream_OFDM_Temp = QAM.Symbol2Bit(x_est_OFDM_Temp);        
            BER_InterferenceCancellation_OFDM(i_SNR,i_rep,i_iterationOFDM) = mean([BinaryDataStream_OFDM_Antenna1;BinaryDataStream_OFDM_Antenna2]~=DetectedBitStream_OFDM_Temp);
        end

        %% Full Block MMSE Equalizer 
        x_est_FullMMSE_OFDM = D_OFDM_MIMO'/(D_OFDM_MIMO*D_OFDM_MIMO'+Pn_time*CorrelationMatrixNoise_OFDM_MIMO)*[y_OFDM_Antenna1(:);y_OFDM_Antenna2(:)];
        x_est_FullMMSE_FBMC = [real(D_FBMC_MIMO);imag(D_FBMC_MIMO)]'/([real(D_FBMC_MIMO);imag(D_FBMC_MIMO)]*[real(D_FBMC_MIMO);imag(D_FBMC_MIMO)]'+Pn_time*CorrelationMatrixNoise_FBMC_RealImag_MIMO)*[real(y_FBMC_Antenna1(:));real(y_FBMC_Antenna2(:));imag(y_FBMC_Antenna1(:));imag(y_FBMC_Antenna2(:))];

        %% Subblock MMSE Equalizer
        x_est_9TapTime_FBMC = FBMC_SubBlockMMSEequalizer2x2MIMO(y_FBMC,D_FBMC_MIMO,'9-tap-Time',Pn_time*CorrelationMatrixNoise_FBMC_RealImag_MIMO,FBMC.Nr.Subcarriers,FBMC.Nr.MCSymbols);    
        x_est_9Tap_FBMC = FBMC_SubBlockMMSEequalizer2x2MIMO(y_FBMC,D_FBMC_MIMO,'9-tap',Pn_time*CorrelationMatrixNoise_FBMC_RealImag_MIMO,FBMC.Nr.Subcarriers,FBMC.Nr.MCSymbols);    
      
        %% Detect Bits (Quantization)
        DetectedBitStream_1TapMMSE_FBMC = PAM.Symbol2Bit(real(x_est_1TapMMSE_FBMC));
        DetectedBitStream_1TapMMSE_OFDM = QAM.Symbol2Bit(x_est_1TapMMSE_OFDM);

        DetectedBitStream_FullMMSE_FBMC = PAM.Symbol2Bit(real(x_est_FullMMSE_FBMC(:)));
        DetectedBitStream_FullMMSE_OFDM = QAM.Symbol2Bit(x_est_FullMMSE_OFDM(:));

        DetectedBitStream_9TapTime_FBMC = PAM.Symbol2Bit(real(x_est_9TapTime_FBMC(:)));
        DetectedBitStream_9Tap_FBMC = PAM.Symbol2Bit(real(x_est_9Tap_FBMC(:)));
       
        %% Calculate Bit Error Ratio
        BER_1TapMMSE_FBMC(i_SNR,i_rep) = mean([BinaryDataStream_FBMC_Antenna1;BinaryDataStream_FBMC_Antenna2]~=DetectedBitStream_1TapMMSE_FBMC);
        BER_1TapMMSE_OFDM(i_SNR,i_rep) = mean([BinaryDataStream_OFDM_Antenna1;BinaryDataStream_OFDM_Antenna2]~=DetectedBitStream_1TapMMSE_OFDM);

        BER_FullMMSE_FBMC(i_SNR,i_rep) = mean([BinaryDataStream_FBMC_Antenna1;BinaryDataStream_FBMC_Antenna2]~=DetectedBitStream_FullMMSE_FBMC);
        BER_FullMMSE_OFDM(i_SNR,i_rep) = mean([BinaryDataStream_OFDM_Antenna1;BinaryDataStream_OFDM_Antenna2]~=DetectedBitStream_FullMMSE_OFDM);

        BER_9TapTimeMMSE_FBMC(i_SNR,i_rep) = mean([BinaryDataStream_FBMC_Antenna1;BinaryDataStream_FBMC_Antenna2]~=DetectedBitStream_9TapTime_FBMC);
        BER_9TapMMSE_FBMC(i_SNR,i_rep) = mean([BinaryDataStream_FBMC_Antenna1;BinaryDataStream_FBMC_Antenna2]~=DetectedBitStream_9Tap_FBMC);

    end
    TimeNeededSoFar = toc;
    disp([int2str(i_rep/NrRepetitions*100) '% Completed! Time Left: ' int2str(TimeNeededSoFar/i_rep*(NrRepetitions-i_rep)/60) 'min, corresponding to approx. '  int2str(TimeNeededSoFar/i_rep*(NrRepetitions-i_rep)/3600) 'hour'])

end
warning('on','MATLAB:nearlySingularMatrix');

%% Simulate Doubly-Flat
disp('Simulate the BER (doubly-flat) ...');
BER_DoublyFlat_MMSE = nan(length(M_SNR_dB),NrRepetitions2);
for i_rep2 = 1: NrRepetitions2
    BinaryDataStream = randi([0 1],log2(QAM.ModulationOrder),2);
    x = QAM.Bit2Symbol(BinaryDataStream);
    H = (randn(2,2)+1j*randn(2,2))/2;    
    y_noNoise = H*x(:);
    for i_SNR2 = 1:length(M_SNR_dB)
        Pn = 10^(-M_SNR_dB(i_SNR2)/10);
        y = y_noNoise + sqrt(Pn/2)*(randn(2,1)+1j*randn(2,1));
        x_est = H'*(H*H'+eye(2)*Pn)^-1*y;
        BER_DoublyFlat_MMSE(i_SNR2,i_rep2) = mean(QAM.Symbol2Bit(x_est)~=BinaryDataStream(:));
    end
end


%% Plot Results
Linewidth1 = 1;
Markersize1 = 6;

figure(44);
semilogy(M_SNR_dB,mean(BER_DoublyFlat_MMSE,2),'black','Linewidth',Linewidth1);
hold on;
semilogy(M_SNR_dB,mean(BER_1TapMMSE_FBMC,2),'black','Linewidth',Linewidth1);
r1 = semilogy(M_SNR_dB,mean(BER_1TapMMSE_OFDM,2),'red --','Linewidth',Linewidth1);

b1=semilogy(M_SNR_dB,mean(BER_9TapTimeMMSE_FBMC,2),'s-','Color',[0 0 1]*1,'Markersize',Markersize1);
b2=semilogy(M_SNR_dB,mean(BER_9TapMMSE_FBMC,2),'-o','Color',[0 0 1]*1,'Markersize',Markersize1);

semilogy(M_SNR_dB,mean(BER_InterferenceCancellation_FBMC(:,:,1),2) ,'-*','Color',[1 0 0]*1,'Markersize',Markersize1);
semilogy(M_SNR_dB,mean(BER_InterferenceCancellation_FBMC(:,:,2),2) ,'-*','Color',[0 0 1]*1,'Markersize',Markersize1);
semilogy(M_SNR_dB,mean(BER_InterferenceCancellation_FBMC(:,:,3),2) ,'-*','Color',[0 1 1]*0.5,'Markersize',Markersize1);
semilogy(M_SNR_dB,mean(BER_InterferenceCancellation_FBMC(:,:,25),2) ,'-*','Color',[1 0 1]*0.8,'Markersize',Markersize1);

p1=semilogy(M_SNR_dB,mean(BER_FullMMSE_FBMC,2),'-x','Color',[0 1 0]*0.8,'Markersize',Markersize1);

bb1 = plot([nan nan],[nan nan],'-* black','Markersize',Markersize1);

legend([r1 b1 b2 bb1 p1],{'CP-OFDM, One-Tap','Suboptimal MMSE, 9-Tap','MMSE Equalizer, 9-Tap','Interference Cancellation','Full Block MMSE'},'Location','SouthWest');

xlabel('Signal-to-Noise Ratio [dB]');
ylabel('Bit Error Ratio');

