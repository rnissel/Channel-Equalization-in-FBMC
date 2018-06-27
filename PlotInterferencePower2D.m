% =========================================================================   
% (c) 2018 Ronald Nissel, ronald.nissel@gmail.com
% =========================================================================      
% This script allows to reproduce Figure 1 and 2 of "FBMC-OQAM in Doubly-
% Selective Channels: A New Perspective on MMSE Equalization", R.Nissel,
% M.Rupp, R.Marsalek, IEEE SPAWC 2017. In particular it shows which
% neighboring time-frequency positions contribute to the interference.

clear; close all;

%% Parameters
FrequencyRange  = 2;                    % Plot range, total range: *2+1
TimeRange       = 4;                    % Plot range, total range: *2+1
Threshold_Upper = -30;                  % 30
Threshold_Lower = -60;                  % 60

% Channel/FBMC Parameters
Velocity_kmh      = 500;                % Velocity in km/h. Note that [mph]*1.6=[kmh] and [m/s]*3.6=[kmh]
PowerDelayProfile = 'VehicularA';       % Channel model, either string or vector: 'Flat', 'AWGN', 'PedestrianA', 'PedestrianB', 'VehicularA', 'VehicularB', 'ExtendedPedestrianA', 'ExtendedPedestrianB', or 'TDL-A_xxns','TDL-B_xxns','TDL-C_xxns' (with xx the RMS delay spread in ns, e.g. 'TDL-A_30ns'), or [1 0 0.2] (Self-defined power delay profile which depends on the sampling rate) 
F                 = 15e3;               % Subcarrier spacing in Hz, 15kHz, same as in LTE
SamplingRate      = F*14*14;            % Sampling rate in Hz. Must be a multiple of the subcarrier spacing. 14 because of the CP in OFDM. F*14*14 because the sampling rate should matche approximately the predefined channel delay taps (Vehicular A)
PrototypeFilter   = 'Hermite';          % Prototype filter for FBMC, either "Hermite" or "PHYDYAS"


%% For Figure 1 
% PrototypeFilter = 'PHYDYAS';


%% FBMC Object
FBMC = Modulation.FBMC(...
    FrequencyRange*2+1,...              % Number subcarriers
    TimeRange*2+1,...                   % Number FBMC symbols
    F,...                               % Subcarrier spacing (Hz)
    SamplingRate,...                    % Sampling rate (Samples/s)
    0,...                               % Intermediate frequency first subcarrier (Hz)
    false,...                           % Transmit real valued signal
    [PrototypeFilter '-OQAM'],...       % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM, 
    4, ...                              % Overlapping factor (also determines oversampling in the frequency domain)
    0, ...                              % Initial phase shift
    true ...                            % Polyphase implementation
    );

%% Channel Model Object
ChannelModel = Channel.FastFading(...
    SamplingRate,...                     % Sampling rate (Samples/s)
    PowerDelayProfile,...                % Power delay profile, either string or vector: 'Flat', 'AWGN', 'PedestrianA', 'PedestrianB', 'VehicularA', 'VehicularB', 'ExtendedPedestrianA', 'ExtendedPedestrianB', or 'TDL-A_xxns','TDL-B_xxns','TDL-C_xxns' (with xx the RMS delay spread in ns, e.g. 'TDL-A_30ns'), or [1 0 0.2] (Self-defined power delay profile which depends on the sampling rate) 
    FBMC.Nr.SamplesTotal,...             % Number of total samples
    Velocity_kmh/3.6*2.5e9/2.998e8,...   % Maximum Doppler shift: Velocity_kmh/3.6*CarrierFrequency/2.998e8  
    'Jakes',...                          % Which Doppler model: 'Jakes', 'Uniform', 'Discrete-Jakes', 'Discrete-Uniform'. For "Discrete-", we assume a discrete Doppler spectrum to improve the simulation time. This only works accuratly if the number of samples and the velocity is sufficiently large                                       
    200, ...                             % Number of paths for the WSSUS process. Only relevant for a 'Jakes' and 'Uniform' Doppler spectrum                                                 
    1,...                                % Number of transmit antennas
    1,...                                % Number of receive antennas
    true ...                             % Gives a warning if the predefined delay taps of the channel do not fit the sampling rate. This is usually not much of a problem if they are approximatly the same.
    );
R_vecH = ChannelModel.GetCorrelationMatrix;

%% Precalculate Transmit and Receive Matrices
G_FBMC = FBMC.GetTXMatrix;
Q_FBMC = (FBMC.GetRXMatrix)';


%% Calculate Interference Power of Surrounding Time-Frequency Positions (reference = middle position)
Position = ceil(FBMC.Nr.Subcarriers/2)+(FBMC.Nr.Subcarriers)*(ceil(FBMC.Nr.MCSymbols/2)-1);
Ey2WithoutDataNoNoise = G_FBMC.'*(kron(sparse(eye(length(Q_FBMC(:,Position)'))),Q_FBMC(:,Position))'*R_vecH*kron(sparse(eye(length(Q_FBMC(:,Position)'))),Q_FBMC(:,Position)))*conj(G_FBMC);
SquarerootEy2 = sqrtm(Ey2WithoutDataNoNoise);
Phase = 1./SquarerootEy2(Position,:).*abs(SquarerootEy2(Position,:));
SquarerootEy2 = SquarerootEy2.*repmat(Phase,[size(SquarerootEy2,1) 1]);
Temp = real(SquarerootEy2)*real(SquarerootEy2)';
P_InterferencePositions = reshape(diag(Temp),FBMC.Nr.Subcarriers,FBMC.Nr.MCSymbols);
P_InterferencePositions = P_InterferencePositions./P_InterferencePositions(FrequencyRange+1,TimeRange+1);

P_InterferencePositions_dB = 10*log10(P_InterferencePositions);


%% Plot Results
ColorMap = 'parula';
colormap(ColorMap);
CurrentColorMap = flipud(colormap);
NewColorMap(:,1) = interp1(linspace(0,100,size(CurrentColorMap,1)),CurrentColorMap(:,1),0:100);
NewColorMap(:,2) = interp1(linspace(0,100,size(CurrentColorMap,1)),CurrentColorMap(:,2),0:100);
NewColorMap(:,3) = interp1(linspace(0,100,size(CurrentColorMap,1)),CurrentColorMap(:,3),0:100);
NewColorMap = [1 1 1;NewColorMap];

P_InterferencePositions_dB(P_InterferencePositions_dB<-60)=-61;
P_InterferencePositions_dB(P_InterferencePositions_dB>-1) = -61;
colormap(NewColorMap);
imagesc(-TimeRange:TimeRange,-FrequencyRange:FrequencyRange,P_InterferencePositions_dB,[Threshold_Lower Threshold_Upper]);

cb = colorbar;
set(cb,'YTick',[-60 -50 -40 -30]);

hold on;
for i_row = (-FrequencyRange+1):(FrequencyRange)
    plot([-.5-TimeRange,TimeRange+0.5],[i_row-.5,i_row-.5],'k-');
end
for i_column = (-TimeRange+1):(TimeRange)
    plot([i_column-.5,i_column-.5],[-.5-FrequencyRange,FrequencyRange+0.5],'k-');
end
plot([-0.5,0.5],[0.5,-0.5],'black');
plot([-0.5,0.5],[-0.5,0.5],'black');
xlabel('Relative Time Position');
ylabel('Relative Frequency Position');


SIR_dB = 10*log10(1./(sum(P_InterferencePositions(:))-1));
disp(['The SIR is ' int2str(SIR_dB) 'dB.']);




