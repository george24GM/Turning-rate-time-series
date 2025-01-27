function ord_pattern_probabilities=lag_patterns_length3(X,d)
% lag = 1
% length p = 3
T=length(X);
matrix_patterns=[2,1,0;0,2,1;1,2,0;2,0,1;1,0,2;0,1,2];
ord_pattern_probabilities=zeros(6,1);
% X_t  X_(t+d)  X_(t+2d)
for t=1:T-2*d
    [~,index]=sort([X(t),X(t+d),X(t+2*d)]);
    index=(flip(index)-1);
    %disp(index)
    ord_number=1;
    found=0;
    while(found==0)
        if (norm(matrix_patterns(ord_number,:)-index)==0)
            ord_pattern_probabilities(ord_number)=ord_pattern_probabilities(ord_number)+1;
            found=1;
        end
       ord_number=ord_number+1;
    end

end
ord_pattern_probabilities=1/(T-(2*d))*ord_pattern_probabilities;
end
