% =========================================================================   
% (c) 2018 Ronald Nissel, ronald.nissel@gmail.com
% =========================================================================      
% This script allows to reproduce Figure 4, 5 and 6 of "FBMC-OQAM in Doubly
% -Selective Channels: A New Perspective on MMSE Equalization", R.Nissel,
% M.Rupp, R.Marsalek, IEEE SPAWC 2017. In particular it simulates the BER
% of OFDM and FBMC in doubly-selective channels and compares different MMSE
% equalization methods as well as a low-complexity interference
% cancellation scheme. 

clear; close all;

%% Simulation Parameters
M_SNR_dB            = [10:2:40];           % SNR in dB
NrRepetitions       = 10;                  % Number of Monte Carlo repetitions (channel realizations) over which we average                  

% Channel
Velocity_kmh        = 500;                 % Velocity in km/h. Note that [mph]*1.6=[kmh] and [m/s]*3.6=[kmh]
PowerDelayProfile   = 'VehicularA';        % Channel model, either string or vector: 'Flat', 'AWGN', 'PedestrianA', 'PedestrianB', 'VehicularA', 'VehicularB', 'ExtendedPedestrianA', 'ExtendedPedestrianB', or 'TDL-A_xxns','TDL-B_xxns','TDL-C_xxns' (with xx the RMS delay spread in ns, e.g. 'TDL-A_30ns'), or [1 0 0.2] (Self-defined power delay profile which depends on the sampling rate) 

% FBMC and OFDM
L                   = 24;                  % Number of subcarriers, 24 to keep the simulation time short and because a higher number does not show significant other behavior. 
F                   = 15e3;                % Subcarrier spacing in Hz, 15kHz, same as in LTE
SamplingRate        = F*14*14;             % Sampling rate in Hz. Must be a multiple of the subcarrier spacing. 14 because of the CP in OFDM. F*14*14 because the sampling rate should matche approximately the predefined channel delay taps (Vehicular A)
QAM_ModulationOrder = 256;                 % Modulation order, 4, 16, 64, 256, 1024, ...
PrototypeFilter     = 'Hermite';           % Prototype filter for FBMC, either "Hermite" or "PHYDYAS"

NrIterationsInterferenceCancellation = 25; % Number of iterations for the simple interference cancellation scheme

%% Parameters in the paper
% NrRepetitions = 2000;

%% For Figure 4
% PrototypeFilter = 'PHYDYAS'


%% Addional Stuff
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
    [PrototypeFilter '-OQAM'],...       % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM, 
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
    1,...                                % Number of transmit antennas
    1,...                                % Number of receive antennas
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

CorrelationMatrixNoise_FBMC_RealImag = [1/2*(real(Q_FBMC')*real(Q_FBMC')'+imag(Q_FBMC')*imag(Q_FBMC')'),...
                            1/2*(real(Q_FBMC')*imag(Q_FBMC')'-imag(Q_FBMC')*real(Q_FBMC')');...
                            -1/2*(real(Q_FBMC')*imag(Q_FBMC')'-imag(Q_FBMC')*real(Q_FBMC')'),...
                            1/2*(real(Q_FBMC')*real(Q_FBMC')'+imag(Q_FBMC')*imag(Q_FBMC')')]; 

                   
%% Preallocate for Parfor
BER_1Tap_FBMC = nan(length(M_SNR_dB),NrRepetitions);
BER_1Tap_OFDM = nan(length(M_SNR_dB),NrRepetitions);

BER_InterferenceCancellation_FBMC = nan(length(M_SNR_dB),NrRepetitions,NrIterationsInterferenceCancellation);
BER_InterferenceCancellation_OFDM = nan(length(M_SNR_dB),NrRepetitions,NrIterationsInterferenceCancellation);

BER_FullMMSE_FBMC = nan(length(M_SNR_dB),NrRepetitions);
BER_FullMMSE_OFDM = nan(length(M_SNR_dB),NrRepetitions);

BER_3TapT_FBMC = nan(length(M_SNR_dB),NrRepetitions);
BER_3TapF_FBMC = nan(length(M_SNR_dB),NrRepetitions);
BER_5Tap_FBMC = nan(length(M_SNR_dB),NrRepetitions);
BER_9Tap_FBMC = nan(length(M_SNR_dB),NrRepetitions);
BER_5TapTBAD_FBMC = nan(length(M_SNR_dB),NrRepetitions);
BER_9TapTBAD_FBMC = nan(length(M_SNR_dB),NrRepetitions); 
BER_13TapTBAD_FBMC = nan(length(M_SNR_dB),NrRepetitions);
BER_21TapTBAD_FBMC = nan(length(M_SNR_dB),NrRepetitions);
BER_13Tap_FBMC = nan(length(M_SNR_dB),NrRepetitions);
BER_21TapTest_FBMC = nan(length(M_SNR_dB),NrRepetitions);

disp('Starting simulation ... ');
tic    
for i_rep = 1:NrRepetitions
    %% Update Channel
    ChannelModel.NewRealization;
    
    %% Binary Data
    BinaryDataStream_FBMC = randi([0 1],FBMC.Nr.Subcarriers*FBMC.Nr.MCSymbols*log2(PAM.ModulationOrder),1);
    BinaryDataStream_OFDM = randi([0 1],OFDM.Nr.Subcarriers*OFDM.Nr.MCSymbols*log2(QAM.ModulationOrder),1);
 
    %% Transmitted Data Symbols (Map bin to symbol)
    x_FBMC = reshape(PAM.Bit2Symbol(BinaryDataStream_FBMC),FBMC.Nr.Subcarriers,FBMC.Nr.MCSymbols);
    x_OFDM = reshape(QAM.Bit2Symbol(BinaryDataStream_OFDM),OFDM.Nr.Subcarriers,OFDM.Nr.MCSymbols);
        
    %% Transmitted FBMC Signal (time domain)
    s_FBMC = G_FBMC*x_FBMC(:); % Same as "FBMC.Modulation(x_FBMX)" which is computationally more efficient. But G_FBMC is consistent with the paper.
    s_OFDM = G_OFDM*x_OFDM(:);
       
    %% Channel
    ConvolutionMatrix = ChannelModel.GetConvolutionMatrix{1};
   
    r_FBMC_noNoise = ConvolutionMatrix*s_FBMC;
    r_OFDM_noNoise = ConvolutionMatrix*s_OFDM;
    
    %% Transmission Matrix
    D_FBMC = Q_FBMC'*ConvolutionMatrix*G_FBMC;
    D_OFDM = Q_OFDM'*ConvolutionMatrix*G_OFDM;

    %% Channel for one tap equalizer
    h_FBMC = diag(D_FBMC);
    h_OFDM = diag(D_OFDM);
    
%     parfor i_SNR = 1:length(M_SNR_dB)      
    for i_SNR = 1:length(M_SNR_dB)
        %% Add Noise
        SNR_dB = M_SNR_dB(i_SNR);
        Pn_time = SamplingRate/(F*L)*10^(-SNR_dB/10);
        noise = sqrt(Pn_time/2)*(randn(size(s_OFDM))+1j*randn(size(s_OFDM)));

        r_FBMC = r_FBMC_noNoise + noise; 
        r_OFDM = r_OFDM_noNoise + noise;

        %% Demodulate FBMC signal
        y_FBMC = reshape(Q_FBMC'*r_FBMC,FBMC.Nr.Subcarriers,FBMC.Nr.MCSymbols); % Same as "FBMC.Demodulation(r_FBMC)" 
        y_OFDM = Q_OFDM'*r_OFDM;

        %% One Tap Equalizer
        x_est_1Tap_FBMC = y_FBMC(:)./h_FBMC;
        x_est_1Tap_OFDM = y_OFDM(:)./h_OFDM;

        %% Simple Interference Cancellation (without sorting etc.)
        x_est_FBMC_Temp = x_est_1Tap_FBMC; % initialize with one tap estimates    
        x_est_OFDM_Temp = x_est_1Tap_OFDM; % initialize with one tap estimates
        for i_iteration = 1:NrIterationsInterferenceCancellation
            y_FBMC_InterferenceCancellation = (y_FBMC(:) - (D_FBMC-diag(h_FBMC))*PAM.SymbolQuantization(real(x_est_FBMC_Temp)));                
            y_OFDM_InterferenceCancellation = (y_OFDM - (D_OFDM-diag(h_OFDM))*QAM.SymbolQuantization(x_est_OFDM_Temp));        

            x_est_FBMC_Temp = y_FBMC_InterferenceCancellation(:)./h_FBMC;
            x_est_OFDM_Temp =  y_OFDM_InterferenceCancellation./h_OFDM;

            DetectedBitStream_FBMC_Temp = PAM.Symbol2Bit(x_est_FBMC_Temp);        
            DetectedBitStream_OFDM_Temp = QAM.Symbol2Bit(x_est_OFDM_Temp);

            BER_InterferenceCancellation_FBMC(i_SNR,i_rep,i_iteration) = mean(BinaryDataStream_FBMC~=DetectedBitStream_FBMC_Temp);
            BER_InterferenceCancellation_OFDM(i_SNR,i_rep,i_iteration) = mean(BinaryDataStream_OFDM~=DetectedBitStream_OFDM_Temp);
        end

        %% Full Block MMSE Equalizer
        x_est_FullMMSE_OFDM = D_OFDM'/(D_OFDM*D_OFDM'+Pn_time*CorrelationMatrixNoise_OFDM)*y_OFDM(:);
        x_est_FullMMSE_FBMC = [real(D_FBMC);imag(D_FBMC)]'/([real(D_FBMC);imag(D_FBMC)]*[real(D_FBMC);imag(D_FBMC)]'+Pn_time*CorrelationMatrixNoise_FBMC_RealImag)*[real(y_FBMC(:));imag(y_FBMC(:))];


        %% Subblock MMSE Equalizer (only for FBMC)
        if strcmp(PrototypeFilter,'Hermite')
            x_est_3TapF_FBMC = FBMC_SubBlockMMSEequalizerHermite(y_FBMC,D_FBMC,'3-tap-frequency',Pn_time*CorrelationMatrixNoise_FBMC_RealImag);    
            x_est_5Tap_FBMC  = FBMC_SubBlockMMSEequalizerHermite(y_FBMC,D_FBMC,'5-tap',Pn_time*CorrelationMatrixNoise_FBMC_RealImag);
            x_est_9Tap_FBMC  = FBMC_SubBlockMMSEequalizerHermite(y_FBMC,D_FBMC,'9-tap',Pn_time*CorrelationMatrixNoise_FBMC_RealImag);
            x_est_13Tap_FBMC = FBMC_SubBlockMMSEequalizerHermite(y_FBMC,D_FBMC,'13-tap',Pn_time*CorrelationMatrixNoise_FBMC_RealImag);

            % Taps only in the time-domain
            x_est_3TapT_FBMC     = FBMC_SubBlockMMSEequalizerHermite(y_FBMC,D_FBMC,'3-tap-time',Pn_time*CorrelationMatrixNoise_FBMC_RealImag);                    
            x_est_5TapTBAD_FBMC  = FBMC_SubBlockMMSEequalizerHermite(y_FBMC,D_FBMC,'5-tap-BAD',Pn_time*CorrelationMatrixNoise_FBMC_RealImag);    
            x_est_9TapTBAD_FBMC  = FBMC_SubBlockMMSEequalizerHermite(y_FBMC,D_FBMC,'9-tap-BAD',Pn_time*CorrelationMatrixNoise_FBMC_RealImag);    
            x_est_13TapTBAD_FBMC = FBMC_SubBlockMMSEequalizerHermite(y_FBMC,D_FBMC,'13-tap-BAD',Pn_time*CorrelationMatrixNoise_FBMC_RealImag);    

        elseif strcmp(PrototypeFilter,'PHYDYAS')      
            x_est_9Tap_FBMC      = FBMC_SubBlockMMSEequalizer(y_FBMC,D_FBMC,'9-tap',Pn_time*CorrelationMatrixNoise_FBMC_RealImag);
            x_est_13Tap_FBMC     = FBMC_SubBlockMMSEequalizer(y_FBMC,D_FBMC,'13-tap',Pn_time*CorrelationMatrixNoise_FBMC_RealImag);
            x_est_21TapTest_FBMC = FBMC_SubBlockMMSEequalizer(y_FBMC,D_FBMC,'21-tap-Test',Pn_time*CorrelationMatrixNoise_FBMC_RealImag);    

            % Taps only in the time-domain
            x_est_3TapT_FBMC     = FBMC_SubBlockMMSEequalizer(y_FBMC,D_FBMC,'3-tap-time',Pn_time*CorrelationMatrixNoise_FBMC_RealImag);            
            x_est_5TapTBAD_FBMC  = FBMC_SubBlockMMSEequalizer(y_FBMC,D_FBMC,'5-tap-BAD',Pn_time*CorrelationMatrixNoise_FBMC_RealImag);    
            x_est_9TapTBAD_FBMC  = FBMC_SubBlockMMSEequalizer(y_FBMC,D_FBMC,'9-tap-BAD',Pn_time*CorrelationMatrixNoise_FBMC_RealImag);    
            x_est_13TapTBAD_FBMC = FBMC_SubBlockMMSEequalizer(y_FBMC,D_FBMC,'13-tap-BAD',Pn_time*CorrelationMatrixNoise_FBMC_RealImag);    
            x_est_21TapTBAD_FBMC = FBMC_SubBlockMMSEequalizer(y_FBMC,D_FBMC,'21-tap-BAD',Pn_time*CorrelationMatrixNoise_FBMC_RealImag);   
        else
            error('Protoype filter must be "Hermite" or "PHYDYAS"');
        end

        %% Detect Bits (Quantization)
        DetectedBitStream_1Tap_FBMC = PAM.Symbol2Bit(real(x_est_1Tap_FBMC(:)));
        DetectedBitStream_1Tap_OFDM = QAM.Symbol2Bit(x_est_1Tap_OFDM(:));

        DetectedBitStream_FullMMSE_FBMC = PAM.Symbol2Bit(real(x_est_FullMMSE_FBMC(:)));
        DetectedBitStream_FullMMSE_OFDM = QAM.Symbol2Bit(x_est_FullMMSE_OFDM(:));

        if strcmp(PrototypeFilter,'Hermite')
            DetectedBitStream_3TapF_FBMC = PAM.Symbol2Bit(real(x_est_3TapF_FBMC(:)));
            DetectedBitStream_5Tap_FBMC  = PAM.Symbol2Bit(real(x_est_5Tap_FBMC(:)));
            DetectedBitStream_9Tap_FBMC  = PAM.Symbol2Bit(real(x_est_9Tap_FBMC(:)));
            DetectedBitStream_13Tap_FBMC  = PAM.Symbol2Bit(real(x_est_13Tap_FBMC(:))); 

            DetectedBitStream_3TapT_FBMC = PAM.Symbol2Bit(real(x_est_3TapT_FBMC(:)));
            DetectedBitStream_5TapTBAD_FBMC  = PAM.Symbol2Bit(real(x_est_5TapTBAD_FBMC(:)));
            DetectedBitStream_9TapTBAD_FBMC  = PAM.Symbol2Bit(real(x_est_9TapTBAD_FBMC(:)));
            DetectedBitStream_13TapTBAD_FBMC  = PAM.Symbol2Bit(real(x_est_13TapTBAD_FBMC(:)));              
        else      
            DetectedBitStream_9Tap_FBMC  = PAM.Symbol2Bit(real(x_est_9Tap_FBMC(:)));
            DetectedBitStream_13Tap_FBMC = PAM.Symbol2Bit(real(x_est_13Tap_FBMC(:)));           
            DetectedBitStream_21TapTest_FBMC  = PAM.Symbol2Bit(real(x_est_21TapTest_FBMC(:)));

            DetectedBitStream_3TapT_FBMC = PAM.Symbol2Bit(real(x_est_3TapT_FBMC(:)));
            DetectedBitStream_5TapTBAD_FBMC  = PAM.Symbol2Bit(real(x_est_5TapTBAD_FBMC(:)));
            DetectedBitStream_9TapTBAD_FBMC  = PAM.Symbol2Bit(real(x_est_9TapTBAD_FBMC(:)));
            DetectedBitStream_13TapTBAD_FBMC  = PAM.Symbol2Bit(real(x_est_13TapTBAD_FBMC(:)));
            DetectedBitStream_21TapTBAD_FBMC  = PAM.Symbol2Bit(real(x_est_21TapTBAD_FBMC(:)));          
        end

        %% Calculate the Bit Error Ratio
        BER_1Tap_FBMC(i_SNR,i_rep) = mean(BinaryDataStream_FBMC~=DetectedBitStream_1Tap_FBMC);
        BER_1Tap_OFDM(i_SNR,i_rep) = mean(BinaryDataStream_OFDM~=DetectedBitStream_1Tap_OFDM);

        BER_FullMMSE_FBMC(i_SNR,i_rep) = mean(BinaryDataStream_FBMC~=DetectedBitStream_FullMMSE_FBMC);
        BER_FullMMSE_OFDM(i_SNR,i_rep) = mean(BinaryDataStream_OFDM~=DetectedBitStream_FullMMSE_OFDM);

        if strcmp(PrototypeFilter,'Hermite')        
            BER_3TapF_FBMC(i_SNR,i_rep) = mean(BinaryDataStream_FBMC~=DetectedBitStream_3TapF_FBMC);
            BER_5Tap_FBMC(i_SNR,i_rep)  = mean(BinaryDataStream_FBMC~=DetectedBitStream_5Tap_FBMC);
            BER_9Tap_FBMC(i_SNR,i_rep)  = mean(BinaryDataStream_FBMC~=DetectedBitStream_9Tap_FBMC);
            BER_13Tap_FBMC(i_SNR,i_rep) = mean(BinaryDataStream_FBMC~=DetectedBitStream_13Tap_FBMC);        

            BER_3TapT_FBMC(i_SNR,i_rep)     = mean(BinaryDataStream_FBMC~=DetectedBitStream_3TapT_FBMC);
            BER_5TapTBAD_FBMC(i_SNR,i_rep)  = mean(BinaryDataStream_FBMC~=DetectedBitStream_5TapTBAD_FBMC);
            BER_9TapTBAD_FBMC(i_SNR,i_rep)  = mean(BinaryDataStream_FBMC~=DetectedBitStream_9TapTBAD_FBMC);
            BER_13TapTBAD_FBMC(i_SNR,i_rep) = mean(BinaryDataStream_FBMC~=DetectedBitStream_13TapTBAD_FBMC);
        else      
            BER_9Tap_FBMC(i_SNR,i_rep)  = mean(BinaryDataStream_FBMC~=DetectedBitStream_9Tap_FBMC);
            BER_13Tap_FBMC(i_SNR,i_rep) = mean(BinaryDataStream_FBMC~=DetectedBitStream_13Tap_FBMC);       
            BER_21TapTest_FBMC(i_SNR,i_rep) = mean(BinaryDataStream_FBMC~=DetectedBitStream_21TapTest_FBMC);

            BER_3TapT_FBMC(i_SNR,i_rep)     = mean(BinaryDataStream_FBMC~=DetectedBitStream_3TapT_FBMC);
            BER_5TapTBAD_FBMC(i_SNR,i_rep)  = mean(BinaryDataStream_FBMC~=DetectedBitStream_5TapTBAD_FBMC);
            BER_9TapTBAD_FBMC(i_SNR,i_rep)  = mean(BinaryDataStream_FBMC~=DetectedBitStream_9TapTBAD_FBMC);
            BER_13TapTBAD_FBMC(i_SNR,i_rep) = mean(BinaryDataStream_FBMC~=DetectedBitStream_13TapTBAD_FBMC);
            BER_21TapTBAD_FBMC(i_SNR,i_rep) = mean(BinaryDataStream_FBMC~=DetectedBitStream_21TapTBAD_FBMC);      
        end    

    end
    TimeNeededSoFar = toc;
    disp([int2str(i_rep/NrRepetitions*100) '% Completed. Time Left: ' int2str(TimeNeededSoFar/i_rep*(NrRepetitions-i_rep)/60) 'min, corresponding to approx. '  int2str(TimeNeededSoFar/i_rep*(NrRepetitions-i_rep)/3600) 'hour'])

end
warning('on','MATLAB:nearlySingularMatrix');


%% Theoretical BEP for a Doubly Flat Rayleigh Channel
M_SNR_dB_morePoints = min(M_SNR_dB):0.5:max(M_SNR_dB);
BitErrorProbability = BitErrorProbabilityDoublyFlatRayleigh(M_SNR_dB_morePoints,PAM.SymbolMapping/sqrt(2),PAM.BitMapping);


%% Plot results
Linewidth1 = 1;
Markersize1 = 6;

if strcmp(PrototypeFilter,'Hermite')   
    figure(5);
    semilogy(M_SNR_dB_morePoints,BitErrorProbability','black','Linewidth',Linewidth1);
    hold on;
    semilogy(M_SNR_dB,nanmean(BER_1Tap_FBMC,2),'black','Linewidth',Linewidth1);
    r1 = semilogy(M_SNR_dB,nanmean(BER_1Tap_OFDM,2),'red --','Linewidth',Linewidth1);

    semilogy(M_SNR_dB,mean(BER_3TapT_FBMC,2),'-s','Color',[1 0 1]*0.8,'Markersize',Markersize1);
    semilogy(M_SNR_dB,mean(BER_5TapTBAD_FBMC,2),'-','Color',[1 1 0]*0.7,'Markersize',Markersize1);
    semilogy(M_SNR_dB(1:2:end),mean(BER_5TapTBAD_FBMC(1:2:end,:),2),'s','Color',[1 1 0]*0.7,'Markersize',Markersize1);
    semilogy(M_SNR_dB,mean(BER_9TapTBAD_FBMC,2),'-','Color',[0 0 1]*1,'Markersize',Markersize1);
    semilogy(M_SNR_dB(1:2:end),mean(BER_9TapTBAD_FBMC(1:2:end,:),2),'s','Color',[0 0 1]*1,'Markersize',Markersize1);
    semilogy(M_SNR_dB,mean(BER_13TapTBAD_FBMC,2),'-','Color',[1 0 0]*1,'Markersize',Markersize1);
    semilogy(M_SNR_dB(1:2:end),mean(BER_13TapTBAD_FBMC(1:2:end,:),2),'s','Color',[1 0 0]*1,'Markersize',Markersize1);

    semilogy(M_SNR_dB,mean(BER_3TapF_FBMC,2),'-','Color',[1 0 1]*0.8,'Markersize',Markersize1);
    semilogy(M_SNR_dB(2:2:end),mean(BER_3TapF_FBMC(2:2:end,:),2),'o','Color',[1 0 1]*0.8,'Markersize',Markersize1);
    semilogy(M_SNR_dB,mean(BER_5Tap_FBMC,2),'-o','Color',[1 1 0]*0.7,'Markersize',Markersize1);
    semilogy(M_SNR_dB,mean(BER_9Tap_FBMC,2),'-o','Color',[0 0 1]*1,'Markersize',Markersize1);
    semilogy(M_SNR_dB,mean(BER_13Tap_FBMC,2),'-o','Color',[1 0 0]*1,'Markersize',Markersize1);
    semilogy(M_SNR_dB,mean(BER_FullMMSE_FBMC,2),'-x','Color',[0 1 0]*0.8,'Markersize',Markersize1);

    b1 = plot([nan nan],[nan nan],'-s black','Markersize',Markersize1);
    b2 = plot([nan nan],[nan nan],'-o black','Markersize',Markersize1);
    b3 = plot([nan nan],[nan nan],'-x','Color',[0 1 0]*0.8,'Markersize',Markersize1);

    xlabel('Signal-to-Noise Ratio [dB]');
    ylabel('Bit Error Ratio');

    legend([r1 b1 b2 b3],{'CP-OFDM, One-Tap','Suboptimal MMSE, n-Tap','MMSE Equalizer, n-Tap','Full Block MMSE'},'Location','SouthWest');

else
    figure(4);
    semilogy(M_SNR_dB_morePoints,BitErrorProbability','black','Linewidth',Linewidth1);
    hold on;
    semilogy(M_SNR_dB,nanmean(BER_1Tap_FBMC,2),'black','Linewidth',Linewidth1);
    r1 = semilogy(M_SNR_dB,nanmean(BER_1Tap_OFDM,2),'red --','Linewidth',Linewidth1);

    semilogy(M_SNR_dB,mean(BER_3TapT_FBMC,2),'Color',[1 0 1]*0.8);
    semilogy(M_SNR_dB,mean(BER_5TapTBAD_FBMC,2),'Color',[1 1 0]*0.7);

    semilogy(M_SNR_dB,mean(BER_9TapTBAD_FBMC,2),'-s','Color',[0 0 1]*1,'Markersize',Markersize1);
    semilogy(M_SNR_dB,mean(BER_13TapTBAD_FBMC,2),'-s','Color',[1 0 0]*1,'Markersize',Markersize1);
    semilogy(M_SNR_dB,mean(BER_21TapTBAD_FBMC,2),'-s','Color',[0 1 1]*0.5,'Markersize',Markersize1);

    semilogy(M_SNR_dB,mean(BER_9Tap_FBMC,2),'-o','Color',[0 0 1]*1,'Markersize',Markersize1);
    semilogy(M_SNR_dB,mean(BER_13Tap_FBMC,2),'-o','Color',[1 0 0]*1,'Markersize',Markersize1);
    semilogy(M_SNR_dB,mean(BER_21TapTest_FBMC,2),'-o','Color',[0 1 1]*0.5,'Markersize',Markersize1);
    semilogy(M_SNR_dB,mean(BER_FullMMSE_FBMC,2),'-x','Color',[0 1 0]*0.8,'Markersize',Markersize1);

    b1 = plot([nan nan],[nan nan],'-s black','Markersize',Markersize1);
    b2 = plot([nan nan],[nan nan],'-o black','Markersize',Markersize1);
    b3 = plot([nan nan],[nan nan],'-x','Color',[0 1 0]*0.8,'Markersize',Markersize1);

    xlabel('Signal-to-Noise Ratio [dB]');
    ylabel('Bit Error Ratio');

    legend([r1 b1 b2 b3],{'CP-OFDM, One-Tap','Suboptimal MMSE, n-Tap','MMSE Equalizer, n-Tap','Full Block MMSE'},'Location','SouthWest');
end


figure(6);
semilogy(M_SNR_dB_morePoints,BitErrorProbability','black','Linewidth',Linewidth1);
hold on;
semilogy(M_SNR_dB,mean(BER_1Tap_FBMC,2),'black','Linewidth',Linewidth1);
r1 = semilogy(M_SNR_dB,mean(BER_1Tap_OFDM,2),'red --','Linewidth',Linewidth1);

g1=semilogy(M_SNR_dB,mean(BER_3TapF_FBMC,2),'-o','Color',[1 1 1]*0.85,'Markersize',Markersize1);
semilogy(M_SNR_dB,mean(BER_5Tap_FBMC,2),'-o','Color',[1 1 1]*0.85,'Markersize',Markersize1);
semilogy(M_SNR_dB,mean(BER_9Tap_FBMC,2),'-o','Color',[1 1 1]*0.85,'Markersize',Markersize1);
semilogy(M_SNR_dB,mean(BER_13Tap_FBMC,2),'-o','Color',[1 1 1]*0.85,'Markersize',Markersize1);

semilogy(M_SNR_dB,mean(BER_InterferenceCancellation_FBMC(:,:,1),2) ,'-*','Color',[1 0 0]*1,'Markersize',Markersize1);
semilogy(M_SNR_dB,mean(BER_InterferenceCancellation_FBMC(:,:,2),2) ,'-*','Color',[0 0 1]*1,'Markersize',Markersize1);
semilogy(M_SNR_dB,mean(BER_InterferenceCancellation_FBMC(:,:,3),2) ,'-*','Color',[0 1 1]*0.5,'Markersize',Markersize1);
semilogy(M_SNR_dB,mean(BER_InterferenceCancellation_FBMC(:,:,25),2) ,'-*','Color',[1 0 1]*0.8,'Markersize',Markersize1);

b1 = plot([nan nan],[nan nan],'-* black','Markersize',Markersize1);

legend([r1 g1 b1],{'CP-OFDM, One-Tap','MMSE Equalizer, 3,5,9,13-Tap','Interference Cancellation'},'Location','SouthWest');
xlabel('Signal-to-Noise Ratio [dB]');
ylabel('Bit Error Ratio');




