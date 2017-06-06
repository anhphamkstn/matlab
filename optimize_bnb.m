function [X,FVAL,EXITFLAG,t,c] = optimize_bnb( Q,q,n,e )
    
    fun = @(x)(1/2)*transpose(x)*Q*x + transpose(q)*x;
    
    eigs_q = eig(Q);
    
    eig_q = eigs_q(1)/2;
    
    ones_n = ones(n,1);
    x0 = ones(n,1)*1/n;
    
    A=[];
    b=[];
    
    Aeq = [ones(1,n);repmat(zeros(1,n),n-1,1)];
    beq = ones(1,n);
    
    xlb = zeros(n);
    xub = ones(n);
    
    g = @(x)(1/2)*transpose(x)*(Q - eig_q*ones_n)*x + transpose(q)*x + (eig_q/2)*transpose(x)*ones(n,1);
    
    z_incumbent=inf; %gia tri thap nhat hien tai
    x_incumbent=inf*ones(n); %nghiem toi uu hien tai
    c = 0; % c la so buoc lap
    t = 0; % t la thoi gian chay
    tic;
    
    
    domain0 = struct('v',ones(n),'x0',x0,'ub_f',fun(x0),'lb_f',fun(x0),'A',A,'b',b); %cau truc mien luu tru gia tri cac canh
   
    [x_temp, z_temp, exitflag]=fmincon(g,x0,A,b,Aeq,beq,xlb,xub);

    if (exitflag == 1) 
        domain0.ub_f = z_temp;domain0.x0 = x_temp;
        domain0.lb_f = fun(x_temp);
        z_incumbent = fun(x_temp);
        x_incumbent = x_temp;
    end
    
    %v là ma tr?n ğ?nh m?i hàng là t?a ğ? m?t ğ?nh
    
    stack = [domain0];
    
    while size(stack) > 0
        c = c+1; 
        %domain_split = struct('v',ones(n),'x',x0,'lamda',[]); %ğõn h?nh s? tách
        split_domain = 1; % ch? s? ğõn h?nh s? tách trong stack
        for i=1:size(stack,1) %ch?n ğõn h?nh s? tách
            if stack(i).ub_f > z_incumbent 
                stack(i) = [];
            elseif stack(i).ub_f <= stack(split_domain).ub_f
                split_domain =i;
            end
        end
        
        %x? l? ğõn h?nh ğ? ch?n
        
        domain_split = stack(split_domain);
        
        stack(split_domain) = [];
        
        %ch?n c?nh dài nh?t
        
        longnest_edge = struct('v1',1,'v2',2,'distant',inf); %c?nh dài nh?t trong ğõn h?nh n?i ğ?nh v1 và v2.
        
        for i=1:n-1
            for j=2:n
                X = [domain_split.v(i,:);domain_split.v(j,:)];
                d = pdist(X,'euclidean');

                if (longnest_edge.distant < d) 
                    longnest_edge.v1 = i;longnest_edge.v2 = j;
                end
            end
        end
                   
        %t?o hai ğõn h?nh m?i
        
        domain1 = domain_split; %ch?a ğ?nh vi
        domain2 = domain_split; %ch?a ğ?nh vj
        
        new_vetex = (domain_split.v(longnest_edge.v1,:)+domain_split.v(longnest_edge.v2,:))/2;
        
        %x? l? domain ch?a vi
        
        domain1.v(longnest_edge.v2,:) = new_vetex;
        domain2.v(longnest_edge.v1,:) = new_vetex;
        
        %tính ub và lb cho t?ng ğõn h?nh m?i chia
        %domain1
        
        x0 = zeros(n);
        for t = 1:n
            x0 = x0 + domain1.v(t);  
        end
        x0 = x0 /n;
         R = zeros(n);%ma tr?n cách vecto c?nh c?a ğõn h?nh
        %phuong trinh sieu phang di qua cac dinh c?a don hinh
        
        sur_v = domain_split.v;
        sur_v(longnest_edge.v1,:) = zeros(1,n);
        sur_v(longnest_edge.v2,:) = new_vetex;
        
       
        
        for i=1:n-1
            temp = sur_v(i,:) - sur_v(i+1,:);
            R(i,:) = temp;
        end
        temp = sur_v(n,:) - sur_v(1,:);
        R(n,:) = temp;
        
        
        surface_eq = @(x)R*x;
        
        [a0,~,~] = fsolve(surface_eq,transpose(ones_n)); 
        
%       g_temp = @(x)(1/2)*transpose(x)*(Q - eig_q*ones_n)*x + transpose(q)*x + (eig_q/2)*transpose(x)*ones(1,n);
        
        domain1.A = [domain1.A;a0];
        domain1.b = [domain1.A;transpose(a0)*new_vetex];
        
        [domain1.x0, domain1.ub_f, ~]=fmincon(g,x0,domain1.A,domain1.b,Aeq,beq,xlb,xub);
        
        domain1.lb_f = fun(domain1.x0);
        
        if (z_incumbent > domain1.lb_f) 
            z_incumbent = domain1.lb_f;
            x_incumbent = domain1.x0;
        end
        
        if (z_incumbent >= domain1.up_f && (domain1.lb_f-domain1.up_f)/domain1.up_f <e)
            stack = [stack domain1];
        end
        
        %domain2
        
        x0 = zeros(n);
        for t = 1:n
            x0 = x0 + domain2.v(t);  
        end
        x0 = x0 /n;
         
        domain2.A = [domain2.A;-a0];
        domain2.b = [domain2.A;transpose(a0)*new_vetex];
        
        [domain2.x0, domain2.ub_f, exitflag]=fmincon(g,x0,domain2.A,domain2.b,Aeq,beq,xlb,xub);
        
        domain2.lb_f = fun(domain2.x0);
        
        if (z_incumbent > domain2.lb_f) 
            z_incumbent = domain2.lb_f;
            x_incumbent = domain2.x0;
        end
        
        if (z_incumbent >= domain2.up_f && (domain1.lb_f-domain1.up_f)/domain1.up_f <e )
            stack = [stack domain2];
        end
                
        %xóa ğõn h?nh ğ? chia
        
    end
    
        
    X = x_incumbent;
    FVAL = x_incumbent;
    EXITFLAG = exitflag;
    t = toc;
    
end

