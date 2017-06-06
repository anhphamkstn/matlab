function [x,fval,exitflag,output] = optimize_bnb_test

    x0 = [ 1/3 ; 2/3];
    H = [3 4;5 6]; 
    f = [-2; -6];
    A = [1 -1];
    b = [0];
    
    fun = @(x)(1/2)*transpose(x)*H*x + transpose(f)*x;
    Aeq = [1 1 ; 0 0];
    beq = [1 ; 0];
    lb = [0;0];
    ub = [1;1];
    
    [x,fval,exitflag,output] = fmincon(fun,x0,-A,b,Aeq,beq,lb,ub);
    
    stacksize=1;
    
    while stacksize > 0
        
        
    
end
