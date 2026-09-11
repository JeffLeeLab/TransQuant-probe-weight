output=[];
for i=1:max(UserData.dots(:,5));
    output(i)=nnz(UserData.dots(:,5)==i);
end;