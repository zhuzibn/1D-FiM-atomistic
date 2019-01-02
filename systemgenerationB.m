Asim_next=zeros(1,natom,'gpuArray');
Asim_previous=zeros(1,natom,'gpuArray');
for ct2=1:natom
    if mark_(ct2)==1%local atom is TM
        if ct2==natom
            Asim_next(ct2)=0;
        else
            if mark_(ct2+1)==1%next atom is TM
                Asim_next(ct2)=A_TMTM;
            else
                Asim_next(ct2)=A_TMRE;
            end
        end
        
        if ct2==1
            Asim_previous(ct2)=0;
        else
            if mark_(ct2-1)==1%previous atom is TM
                Asim_previous(ct2)=A_TMTM;
            else
                Asim_previous(ct2)=A_TMRE;
            end
        end
    else%local atom is RE
        if ct2==natom
            Asim_next(ct2)=0;
        else
            if mark_(ct2+1)==1%next atom is TM
                Asim_next(ct2)=A_TMRE;
            else
                Asim_next(ct2)=A_RERE;
            end
        end
        
        if ct2==1
            Asim_previous(ct2)=0;
        else
            if mark_(ct2-1)==1%previous atom is TM
                Asim_previous(ct2)=A_TMRE;
            else
                Asim_previous(ct2)=A_RERE;
            end
        end
    end
end
clear ct2