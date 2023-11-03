function Y= extract_function(k,slope,Mx,My,p1)
    x=k;
    for i =1:length(p1)
        A(i)=p1(i)*(x^(length(p1)-i));
    end
    Y=sum(A)-(slope*(x-Mx)+My);

end