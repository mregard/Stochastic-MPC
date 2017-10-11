function U = generateBinarySignal(lims,alpha,N)
U=zeros(1,N);
U(1)=lims(1);
for i=2:N
    p = rand;
    if p<alpha && U(i-1)==lims(1)
        U(i)=lims(2);
    elseif p<alpha
        U(i)=lims(1);
    else
        U(i)=U(i-1);
    end
end


end

