function Y= raffinate_function(k,slope,Mx,My,p2)
    x=k;
    for i =1:length(p2)
        A(i)=p2(i)*(x^(length(p2)-i));
    end
    Y=sum(A)-(slope*(x-Mx)+My);
    

end