function [sgates,errors]=GMR2(signal,numsurrogates,thresh)
%Gradual Multifractal Reconstruction in two dimensions.
%The base algorithm (the IAAWT) was developed in:
%
%Keylock, C.J. 2017. Multifractal surrogate-data generation algorithm that 
%preserves pointwise Hölder regularity structure, with initial applications
%to turbulence, Physical Review E 95, 032123, https://doi.org/10.1103/PhysRevE.95.032123.
%
%GMR was introduced in:
%Keylock, C.J. 2018. Gradual multifractal reconstruction of time-series: 
%Formulation of the method and an application to the coupling between stock
%market indices and their Hölder exponents, Physica D 368, 1-9, https://doi.org/10.1016/j.physd.2017.11.011
%
%An example in two dimensions appeared in
%Keylock, C.J. 2019. Hypothesis testing for nonlinear phenomena in the 
%geosciences using synthetic, surrogate data, Earth and Space Science 6, doi: 10.1029/2018EA000435
%
%And was first used properly in a geomorphological application:
%Keylock, C.J., Singh, A., Passalacqua, P., and Foufoula-Georgiou, E. 2021.
%``Dissecting'' Landscapes with Holder Exponents to Reconcile Process and Form

%The code here is adopted to employ the new functions that exist in the
%Matlab Wavelet Toolbox. Our original code used the original toolbox
%developed by Nick Kingsbury who invented the DTCWT transform:
%
%Kingsbury, N. 2001. Complex wavelets for shift invariant analysis and
% filtering of signals, Appl. Comput. Harmon. Anal. 10, 234-253.

%The user inputs a 2-D signal, which is square with sides 2^N where N is an
%integer. The user states the number of surrogates to generate and the
%threhold (the GMR parameter, which defaults to 0 - the IAAWT algorithm). 

%The code returns a structure containing the surrogates 
%as well as a convergence error measure.

sizer=size(signal);

numlevels=floor(log(min(sizer))/log(2));
if rem(log(length(signal))/log(2),1)~=0
    signal=signal(1:2^numlevels,1:2^numlevels);
end
sizer=size(signal);

no_arg=nargin;
if no_arg<3
    thresh=0;
end
if no_arg<2
   numsurrogates=1;
end

sortval=sort(reshape(signal,sizer(1)*sizer(2),1));
stdval=std(sortval);

for loop1=1:numsurrogates
    sgates{loop1}=zeros(sizer);
end

%Do the transform and work out the scaled energies
%wavelet decomp
%[Yl,Yh] = dtwavexfm2(signal,numlevels,'near_sym_b','qshift_b');
DT = dddtree2('cplxdt',signal,numlevels-3,'dtf3');
%Convert into complex form
for loop1=1:length(DT.cfs)-1
    Yh{loop1}=squeeze(DT.cfs{loop1}(:,:,:,:,1))+i*squeeze(DT.cfs{loop1}(:,:,:,:,2));
end
Yh{length(DT.cfs)}=squeeze(DT.cfs{length(DT.cfs)}(:,:,:,1))+i*squeeze(DT.cfs{length(DT.cfs)}(:,:,:,2));
    
startval=1;
Amplitudes=zeros(2^((numlevels-1)*2)*6,1);
for loop1=1:numlevels-3;
    ampYh{loop1,1}=abs(Yh{loop1});
    sizeval{loop1}=size(ampYh{loop1});
    scalesAmps{loop1}=(abs(Yh{loop1}).^2)./(2^loop1);
    endval=startval+(sizeval{loop1}(1)*sizeval{loop1}(2)*6)-1;
    Amplitudes(startval:endval)=reshape(scalesAmps{loop1},sizeval{loop1}(1)*sizeval{loop1}(2)*6,1);
    startval=endval+1;
end
loop1=numlevels-2;
ampYh{loop1,1}=abs(Yh{loop1});
sizeval{loop1}=size(ampYh{loop1});
scalesAmps{loop1}=(abs(Yh{loop1}).^2)./(2^loop1);
endval=startval+(sizeval{loop1}(1)*sizeval{loop1}(2)*2)-1;
Amplitudes(startval:endval)=reshape(scalesAmps{loop1},sizeval{loop1}(1)*sizeval{loop1}(2)*2,1);
startval=endval+1;

%What do we fix? 
[sort_wsq,wsq_pos]=sort(Amplitudes,'descend');
summed=cumsum(sort_wsq);
toten=summed(length(summed));
summed=summed/toten;
threshpos=find(summed>thresh,1,'first');
%The threshold applies to the values in scalesAmps
threshold=sort_wsq(threshpos);

%Write a binary look-up where 1=fixed

for loop1=1:numlevels-3
    LookUp{loop1,1}(1:sizeval{loop1}(1),1:sizeval{loop1}(2),1:6)=0;
    temp=find(scalesAmps{loop1}>threshold);
    LookUp{loop1,1}(temp)=1;
end
loop1=numlevels-2;
LookUp{loop1,1}(1:sizeval{loop1}(1),1:sizeval{loop1}(2),1:2)=0;
temp=find(scalesAmps{loop1}>threshold);
LookUp{loop1,1}(temp)=1;

%Now have the architecture to do the constrained randomisation

for surrloop=1:numsurrogates
    disp(strcat('Making surrogate number_',num2str(surrloop)))
    
    %make a random dataset and take its imag and phases
    [dummy,shuffind]=sort(rand(size(sortval)));
    linearz(shuffind,1) = sortval;
    z=reshape(linearz,sizer(1),sizer(2));
    
    %[Zl,Zh] = dtwavexfm2(z,maxlevels,'near_sym_b','qshift_b');
    ZT = dddtree2('cplxdt',z,numlevels-3,'dtf3');
    for loopz=1:length(ZT.cfs)-1
        Zh{loopz}=squeeze(ZT.cfs{loopz}(:,:,:,:,1))+i.*squeeze(ZT.cfs{loopz}(:,:,:,:,2));
    end
    Zh{length(ZT.cfs)}=squeeze(ZT.cfs{length(ZT.cfs)}(:,:,:,1))+i.*squeeze(ZT.cfs{length(ZT.cfs)}(:,:,:,2));
    
    for loop1=1:numlevels-2
        newphase{loop1}=angle(Zh{loop1});
    end
    
    accerror=.0001;
    amperror(1)=100;
    waverror(1)=100;   
    counter=1;
    newZh=Yh;
    while (amperror(counter) > accerror) && (waverror(counter) > accerror)
        %wavelet construction
        oldz=z;
        for loop1=1:numlevels-2
            %Find the points where LookUp=0 and randomise these only
            temp=find(LookUp{loop1}==0);
            if ~isempty(temp)
                newZh{loop1}(temp)=ampYh{loop1}(temp).*exp(i.*newphase{loop1}(temp));
            end
        end
        
        %Convert back out from complex form
        for loopz=1:length(DT.cfs)-1
            DT.cfs{loopz}(:,:,:,:,1)=real(newZh{loopz});
            DT.cfs{loopz}(:,:,:,:,2)=imag(newZh{loopz});
        end
        DT.cfs{length(DT.cfs)}(:,:,:,1)=real(newZh{length(DT.cfs)});
        DT.cfs{length(DT.cfs)}(:,:,:,2)=imag(newZh{length(DT.cfs)});

        z=idddtree2(DT);
        %z=dtwaveifm2(Yl,newZh,'near_sym_b','qshift_b');
        wavdiff=mean(mean(abs(real(z)-real(oldz))));
        waverror(counter+1) = wavdiff/stdval;
        
        %impose original values       
        oldz=z;
        lineardata=reshape(z,sizer(1)*sizer(2),1);
        [dummy,shuffind]=sort(real(lineardata));
        linearz(shuffind)=sortval;
        z=reshape(linearz,sizer(1),sizer(2));
        ampdiff=mean(mean(abs(real(z)-real(oldz))));
        amperror(counter+1) = ampdiff/stdval;
        
        %Wavelet step
        ZT = dddtree2('cplxdt',z,numlevels-3,'dtf3');
        for loopz=1:length(ZT.cfs)-1
            Zh{loopz}=squeeze(ZT.cfs{loopz}(:,:,:,:,1))+i.*squeeze(ZT.cfs{loopz}(:,:,:,:,2));
        end
        Zh{length(ZT.cfs)}=squeeze(ZT.cfs{length(ZT.cfs)}(:,:,:,1))+i.*squeeze(ZT.cfs{length(ZT.cfs)}(:,:,:,2));

        %[nZl,nZh] = dtwavexfm2(z,maxlevels,'near_sym_b','qshift_b');
        %get phases and imag
        for loop1=1:length(ZT.cfs)
            newphase{loop1}=angle(Zh{loop1});
        end
        
        toterror=amperror(counter+1)+waverror(counter+1);
        oldtoterr=amperror(counter)+waverror(counter);
        if abs((oldtoterr-toterror)/toterror) < (accerror/10);
            amperror(counter+1)=-1;
        end
        counter=counter+1;
    end
    
    clear amperror waverror
    sgates{surrloop}=z;
    errors(surrloop,1)=toterror;
end