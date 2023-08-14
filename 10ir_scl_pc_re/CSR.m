function regout=CSR(regin)
n=length(regin);
for i=1:n
    regout=[regin(2:n),regin(1)];
end
end