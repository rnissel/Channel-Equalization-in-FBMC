function x_est = FBMC_SubBlockMMSEequalizer2x2MIMO(...
    y, ...       % Received samples
    D, ...       % Transmission matrix: y = D *x + n;
    Method,...   % Method, i.e., how many taps: '3-tap-time','3-tap-frequency', '5-tap', '9-tap', '19-tap'
    Rn,...       % Noise correlation matrix
    L,...        % Number of subcarriers
    K,...        % Number of time symbols
    Rx ...       % Data symbol correlation matrix
    )  

if not(exist('Rx','var'))
    UncorrelatedInput = true;
else
    Rx = sparse(Rx);
end

x_est = nan(L*K*2,1);
for i_lk = 1:L*K
    SubblockMatrix = zeros(L,K);    
    SubblockMatrix(i_lk) = 1;    
    switch Method
         case '9-tap-Time'
            [a,b] = find(SubblockMatrix);
            PositionTime = (b-4):(b+4);
            PositionTime(PositionTime<1)=1;
            PositionTime(PositionTime>K)=K;       
            SubblockMatrix(a,PositionTime)=-1;
            SubblockMatrix(i_lk)=1; 
        case '9-tap'
            [a,b] = find(SubblockMatrix);
            PositionFrequency = (a-1):(a+1);        
            PositionFrequency(PositionFrequency<1)=1;
            PositionFrequency(PositionFrequency>L)=L;                
            PositionTime = (b-1):(b+1);
            PositionTime(PositionTime<1)=1;
            PositionTime(PositionTime>K)=K;    
            SubblockMatrix(PositionFrequency,PositionTime)=-1;   
            SubblockMatrix(i_lk)=1;               
        otherwise
            error('Method not supported');
    end
    
    SubblockMatrix = [SubblockMatrix(:);SubblockMatrix(:)];
        
    yS = y(SubblockMatrix~=0);
    DS = D(SubblockMatrix~=0,SubblockMatrix~=0);
    DR = D(SubblockMatrix~=0,SubblockMatrix==0);        
    RnS = Rn([SubblockMatrix(:)~=0;SubblockMatrix(:)~=0],[SubblockMatrix(:)~=0;SubblockMatrix(:)~=0]);
    
    if UncorrelatedInput    
        DSlk = D(SubblockMatrix~=0,SubblockMatrix==1);            
        x_est(SubblockMatrix==1) = [real(DSlk);imag(DSlk)]'/([real(DS);imag(DS)]*[real(DS);imag(DS)]'+[real(DR);imag(DR)]*[real(DR);imag(DR)]'+RnS)*[real(yS);imag(yS)];         
    else
%         RxS = Rx(SubblockMatrix~=0,SubblockMatrix~=0);
%         asdf = Rx(not(SubblockMatrix),not(SubblockMatrix));
%         RxSR = Rx(SubblockMatrix~=0,SubblockMatrix==0);
%         RxRS = Rx(SubblockMatrix==0,SubblockMatrix~=0);
%     
%         x_est_block = RxS'*[real(DS);imag(DS)]'/([real(DS);imag(DS)]*RxS*[real(DS);imag(DS)]'+[real(DR);imag(DR)]*asdf*[real(DR);imag(DR)]'+[real(DS);imag(DS)]*RxSR*[real(DR);imag(DR)]'+[real(DR);imag(DR)]*RxRS*[real(DS);imag(DS)]'+RnS)*[real(yS);imag(yS)];   
%         x_est(i_lk) = x_est_block(SubblockMatrix(SubblockMatrix~=0)==1);
            error('Not implemented');
    end
   
    
%     imagesc(SubblockMatrix);
%     pause(0.03)
end
        



