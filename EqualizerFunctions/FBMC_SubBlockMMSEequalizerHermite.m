function x_est = FBMC_SubBlockMMSEequalizer(...
    y, ...       % Received samples
    D, ...       % Transmission matrix: y = D *x + n;
    Method,...   % Method, i.e., how many taps: '3-tap-time','3-tap-frequency', '5-tap', '9-tap', '19-tap'
    Rn,...       % Noise correlation matrix
    Rx ...       % Data symbol correlation matrix
    )  

if not(exist('Rx','var'))
    UncorrelatedInput = true;
else
    Rx = sparse(Rx);
end

L = size(y,1);  % Number of "subcarriers"
K = size(y,2);  % Number of "time-symbols"

x_est = nan(L*K,1);
for i_lk = 1:L*K
    SubblockMatrix = zeros(L,K);    
    SubblockMatrix(i_lk) = 1;    
    switch Method
         case '5-tap-BAD'
            [a,b] = find(SubblockMatrix);
            PositionTime = (b-2):(b+2);
            PositionTime(PositionTime<1)=1;
            PositionTime(PositionTime>K)=K;       
            SubblockMatrix(a,PositionTime)=-1;
            SubblockMatrix(i_lk)=1;
         case '9-tap-BAD'
            [a,b] = find(SubblockMatrix);
            PositionTime = (b-4):(b+4);
            PositionTime(PositionTime<1)=1;
            PositionTime(PositionTime>K)=K;       
            SubblockMatrix(a,PositionTime)=-1;
            SubblockMatrix(i_lk)=1;
         case '13-tap-BAD'
            [a,b] = find(SubblockMatrix);
            PositionTime = (b-6):(b+6);
            PositionTime(PositionTime<1)=1;
            PositionTime(PositionTime>K)=K;       
            SubblockMatrix(a,PositionTime)=-1;
            SubblockMatrix(i_lk)=1;                
         case '15-tap-BAD'
            [a,b] = find(SubblockMatrix);
            PositionTime = (b-7):(b+7);
            PositionTime(PositionTime<1)=1;
            PositionTime(PositionTime>K)=K;       
            SubblockMatrix(a,PositionTime)=-1;
            SubblockMatrix(i_lk)=1;   
         case '21-tap-BAD'
            [a,b] = find(SubblockMatrix);
            PositionTime = (b-10):(b+10);
            PositionTime(PositionTime<1)=1;
            PositionTime(PositionTime>K)=K;       
            SubblockMatrix(a,PositionTime)=-1;
            SubblockMatrix(i_lk)=1;              
         case '25-tap-BAD'
            [a,b] = find(SubblockMatrix);
            PositionTime = (b-12):(b+12);
            PositionTime(PositionTime<1)=1;
            PositionTime(PositionTime>K)=K;       
            SubblockMatrix(a,PositionTime)=-1;
            SubblockMatrix(i_lk)=1;               
        case '3-tap-time'
            [a,b] = find(SubblockMatrix);
            PositionTime = (b-1):(b+1);
            PositionTime(PositionTime<1)=1;
            PositionTime(PositionTime>K)=K;       
            SubblockMatrix(a,PositionTime)=-1;
            SubblockMatrix(i_lk)=1;
        case '3-tap-frequency'
            [a,b] = find(SubblockMatrix);
            PositionFrequency = (a-1):(a+1);        
            PositionFrequency(PositionFrequency<1)=1;
            PositionFrequency(PositionFrequency>L)=L;       
            SubblockMatrix(PositionFrequency,b)=-1;
            SubblockMatrix(i_lk)=1;   
        case '5-tap-frequency'
            [a,b] = find(SubblockMatrix);
            PositionFrequency = (a-2):(a+2);        
            PositionFrequency(PositionFrequency<1)=1;
            PositionFrequency(PositionFrequency>L)=L;       
            SubblockMatrix(PositionFrequency,b)=-1;
            SubblockMatrix(i_lk)=1;             
        case '5-tap'
            [a,b] = find(SubblockMatrix);
            PositionFrequency = (a-1):(a+1);        
            PositionFrequency(PositionFrequency<1)=1;
            PositionFrequency(PositionFrequency>L)=L;                
            PositionTime = (b-1):(b+1);
            PositionTime(PositionTime<1)=1;
            PositionTime(PositionTime>K)=K;    
            SubblockMatrix(PositionFrequency,b)=-1;
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
        case '13-tap'
            [a,b] = find(SubblockMatrix);
            PositionFrequency = (a-1):(a+1);        
            PositionFrequency(PositionFrequency<1)=1;
            PositionFrequency(PositionFrequency>L)=L;                
            PositionTime = (b-1):(b+1);
            PositionTime(PositionTime<1)=1;
            PositionTime(PositionTime>K)=K;          
            SubblockMatrix(PositionFrequency,PositionTime)=-1;   
            SubblockMatrix(PositionFrequency,PositionTime)=-1; 
            if b>2
                SubblockMatrix(a,b-2)=-1;
            end
            if b<=K-2
                SubblockMatrix(a,b+2)=-1;               
            end     
            if a>2
                SubblockMatrix(a-2,b)=-1;
            end
            if a<=L-2
                SubblockMatrix(a+2,b)=-1;               
            end   
            SubblockMatrix(i_lk)=1;              
        case '15-tap'
            [a,b] = find(SubblockMatrix);
            PositionFrequency = (a-1):(a+1);        
            PositionFrequency(PositionFrequency<1)=1;
            PositionFrequency(PositionFrequency>L)=L;                
            PositionTime = (b-2):(b+2);
            PositionTime(PositionTime<1)=1;
            PositionTime(PositionTime>K)=K;          
            SubblockMatrix(PositionFrequency,PositionTime)=-1;   
            SubblockMatrix(PositionFrequency,PositionTime)=-1; 
            if b>2
                SubblockMatrix(a,b-2)=0;
            end
            if b<=K-2
                SubblockMatrix(a,b+2)=0;               
            end     
            if a>2
                SubblockMatrix(a-2,b)=-1;
            end
            if a<=L-2
                SubblockMatrix(a+2,b)=-1;               
            end   
            SubblockMatrix(i_lk)=1; 
        case '21-tap'
            [a,b] = find(SubblockMatrix);
            PositionFrequency = (a-2):(a+2);        
            PositionFrequency(PositionFrequency<1)=1;
            PositionFrequency(PositionFrequency>L)=L;                
            PositionTime = (b-2):(b+2);
            PositionTime(PositionTime<1)=1;
            PositionTime(PositionTime>K)=K;          
            SubblockMatrix(PositionFrequency,PositionTime)=-1;   
            SubblockMatrix(PositionFrequency,PositionTime)=-1;   
            if b>2
                if a>2
                    SubblockMatrix(a-2,b-2)=0;
                end
                if a<L-2
                    SubblockMatrix(a+2,b-2)=0;                    
                end
            end  
            if b<=K-2
                if a>2
                    SubblockMatrix(a-2,b+2)=0;
                end
                if a<L-2
                    SubblockMatrix(a+2,b+2)=0;                    
                end
            end              
            SubblockMatrix(i_lk)=1;             
        case '25-tap'
            [a,b] = find(SubblockMatrix);
            PositionFrequency = (a-2):(a+2);        
            PositionFrequency(PositionFrequency<1)=1;
            PositionFrequency(PositionFrequency>L)=L;                
            PositionTime = (b-2):(b+2);
            PositionTime(PositionTime<1)=1;
            PositionTime(PositionTime>K)=K;    
            SubblockMatrix(PositionFrequency,PositionTime)=-1;   
            SubblockMatrix(i_lk)=1;  
        case '13-tap-Sym'
            [a,b] = find(SubblockMatrix);
            PositionFrequency = (a-1):(a+1);        
            PositionFrequency(PositionFrequency<1)=1;
            PositionFrequency(PositionFrequency>L)=L;                
            PositionTime = (b-1):(b+1);
            PositionTime(PositionTime<1)=1;
            PositionTime(PositionTime>K)=K;          
            SubblockMatrix(PositionFrequency,PositionTime)=-1;   
            SubblockMatrix(PositionFrequency,PositionTime)=-1; 
            if b>2
               SubblockMatrix(a,b-2)=-1;  
            end
            if b<=K-2
               SubblockMatrix(a,b+2)=-1;
            end
            if a>2
               SubblockMatrix(a-2,b)=-1;  
            end
            if a<=L-2
               SubblockMatrix(a+2,b)=-1;
            end            
            SubblockMatrix(i_lk)=1;  
        case '25-tap-Sym'
            [a,b] = find(SubblockMatrix);
            PositionFrequency = (a-2):(a+2);        
            PositionFrequency(PositionFrequency<1)=1;
            PositionFrequency(PositionFrequency>L)=L;                
            PositionTime = (b-2):(b+2);
            PositionTime(PositionTime<1)=1;
            PositionTime(PositionTime>K)=K;          
            SubblockMatrix(PositionFrequency,PositionTime)=-1;   
            SubblockMatrix(PositionFrequency,PositionTime)=-1;   
            SubblockMatrix(i_lk)=1;        
        case '35-tap-Test'
            [a,b] = find(SubblockMatrix);
            PositionFrequency = (a-2):(a+2);        
            PositionFrequency(PositionFrequency<1)=1;
            PositionFrequency(PositionFrequency>L)=L;                
            PositionTime = (b-3):(b+3);
            PositionTime(PositionTime<1)=1;
            PositionTime(PositionTime>K)=K;          
            SubblockMatrix(PositionFrequency,PositionTime)=-1;   
            SubblockMatrix(PositionFrequency,PositionTime)=-1;   
            SubblockMatrix(i_lk)=1;   
        case '21-tap-Test'
            [a,b] = find(SubblockMatrix);
            PositionFrequency = (a-1):(a+1);        
            PositionFrequency(PositionFrequency<1)=1;
            PositionFrequency(PositionFrequency>L)=L;                
            PositionTime = (b-2):(b+2);
            PositionTime(PositionTime<1)=1;
            PositionTime(PositionTime>K)=K;          
            SubblockMatrix(PositionFrequency,PositionTime)=-1;   
            SubblockMatrix(PositionFrequency,PositionTime)=-1; 
            
            PositionFrequency = (a-2):(a+2);        
            PositionFrequency(PositionFrequency<1)=1;
            PositionFrequency(PositionFrequency>L)=L;                
            PositionTime = (b-1):(b+1);
            PositionTime(PositionTime<1)=1;
            PositionTime(PositionTime>K)=K;          
            SubblockMatrix(PositionFrequency,PositionTime)=-1;   
            SubblockMatrix(PositionFrequency,PositionTime)=-1;             
            
            if b>3
               SubblockMatrix(a,b-3)=-1;  
            end
            if b<=K-3
               SubblockMatrix(a,b+3)=-1;
            end
            if b>2
                SubblockMatrix(a,b-2)=0;
            end
            if b<=K-2
                SubblockMatrix(a,b+2)=0;               
            end            
            SubblockMatrix(i_lk)=1;  
        case '55-tap'
            [a,b] = find(SubblockMatrix);
            PositionFrequency = (a-2):(a+2);        
            PositionFrequency(PositionFrequency<1)=1;
            PositionFrequency(PositionFrequency>L)=L;                
            PositionTime = (b-5):(b+5);
            PositionTime(PositionTime<1)=1;
            PositionTime(PositionTime>K)=K;          
            SubblockMatrix(PositionFrequency,PositionTime)=-1;   
            SubblockMatrix(i_lk)=1;              
        otherwise
            error('Method not supported');
    end
        
    yS = y(SubblockMatrix~=0);
    DS = D(SubblockMatrix~=0,SubblockMatrix~=0);
    DR = D(SubblockMatrix~=0,SubblockMatrix==0);        
    RnS = Rn([SubblockMatrix(:)~=0;SubblockMatrix(:)~=0],[SubblockMatrix(:)~=0;SubblockMatrix(:)~=0]);
    
    if UncorrelatedInput    
        DSlk = D(SubblockMatrix~=0,SubblockMatrix==1);            
        x_est(i_lk) = [real(DSlk);imag(DSlk)]'/([real(DS);imag(DS)]*[real(DS);imag(DS)]'+[real(DR);imag(DR)]*[real(DR);imag(DR)]'+RnS)*[real(yS);imag(yS)];         
    else
        RxS = Rx(SubblockMatrix~=0,SubblockMatrix~=0);
        asdf = Rx(not(SubblockMatrix),not(SubblockMatrix));
        RxSR = Rx(SubblockMatrix~=0,SubblockMatrix==0);
        RxRS = Rx(SubblockMatrix==0,SubblockMatrix~=0);
    
        x_est_block = RxS'*[real(DS);imag(DS)]'/([real(DS);imag(DS)]*RxS*[real(DS);imag(DS)]'+[real(DR);imag(DR)]*asdf*[real(DR);imag(DR)]'+[real(DS);imag(DS)]*RxSR*[real(DR);imag(DR)]'+[real(DR);imag(DR)]*RxRS*[real(DS);imag(DS)]'+RnS)*[real(yS);imag(yS)];   
        x_est(i_lk) = x_est_block(SubblockMatrix(SubblockMatrix~=0)==1);
    end
   
    
%     imagesc(SubblockMatrix);
%     pause(0.03)
end
        



