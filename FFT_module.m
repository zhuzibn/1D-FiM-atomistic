%% FFT for time domain => freq domain, & locate fundamental freq.
% input:
% 1. nT: number of sampling points
% 2. runTime: time duration, [ns]
% 3. y: time domain information
% 4. rminit:remove percentage of data from the beginning
% 5. rmlast:remove percentage of data from the beginning
% 6. plotfft, plot FFT spectrum (1) or not(0)
function centfreq=FFT_module(nT0,runTime0,y,rminit,rmlast,plotfft)
nT=nT0;
runTime=runTime0;
%plottotal=0; %plot total spectrum or partial

rm_amount_head=rminit; %remove first xx data
nT=ceil(nT-nT0*rm_amount_head);
runTime=runTime-runTime0*rm_amount_head;
y=y(:,(end-nT):end);

rm_amount_tail=rmlast;%remove last 1/8 data
nT=ceil(nT-nT0*rm_amount_tail);
runTime=runTime-runTime0*rm_amount_tail;
y=y(:,1:nT);

L=nT;
Fs=L/(runTime/1e9); %sampling frequency
%% total frequency domain
freq_1=(0:ceil(L/2))*Fs/L;     % frequency vector

tmp1=2*abs(fft(y)/L);%fourier transform of mx, double side spectrum
fftx_1=tmp1(5:ceil(L/2)+1); % single side spectrum

if plotfft
    figure;
    plot(freq_1(5:end)*1e-9,fftx_1)
    xlabel('frequency(GHz)');ylabel('mV');
end

[~,col]=max(fftx_1); % locate fundamental frequency
centfreq=freq_1(col)*1e-9;
if (0)%debug use, visualization
figure('name','FFT-mx')
plot(freq_1(ceil(col/2):col*4)*1e-9,fftx_1(ceil(col/2):col*4))
xlabel('frequency(GHz)');ylabel('mV');
end

if (0) %to do, PSD for y
    tmp1=fft(y);
    tmp1=tmp1(1:ceil(L/2)+1);
    psdd=(1/(L^2))*abs(tmp1).^2;
    psdd(2:end-1)=2*psdd(2:end-1);
    
    %xdft = xdft(1:N/2+1);
    %psdx = (1/(2*pi*N)) * abs(xdft).^2;
    plot(freq_1(2:end-1),10*log10(psdd(2:end-1)),'r')
end
end







