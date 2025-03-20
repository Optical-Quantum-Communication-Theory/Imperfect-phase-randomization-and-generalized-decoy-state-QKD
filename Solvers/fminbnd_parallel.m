%iterative linear search of a function f from xmin to xmax
%can be parallelized with parfor
%usage: stepSize = fminbnd_parallel(@(t)primalf(rho+t*deltaRho,keyMap,krausOperators),0,1,10,3);

function [xopt,fopt] = fminbnd_parallel(f,xmin,xmax,steps,depth)
    N = steps; %number of sample points
    stepsize = (xmax-xmin)/(N-1); %size of each step
    
    %first iteration
    results = zeros(1,N);
    %fprintf('searching %f to %f\n',xmin,xmax)
    parfor i=1:N
        xt = xmin + stepsize*(i-1);
        results(i) = f(xt);
    end
    [fopt,iopt] = min(results);
    xopt=xmin+stepsize*(iopt-1);
    
    depth = depth - 1;
    
    while(depth>0)
        %further iteration (optional)
        %perform another search within one block to the left/right of optimal sample point
        imin = max(1,iopt-1);
        imax = min(N,iopt+1);
        xmax = xmin+stepsize*(imax-1);
        xmin = xmin+stepsize*(imin-1);
        stepsize = (xmax-xmin)/(N-1);
        results = zeros(1,N);
        
        %fprintf('searching %f to %f\n',xmin,xmax)
        parfor i=1:N
            xt = xmin + stepsize*(i-1);
            results(i) = f(xt);
        end
        [fopt,iopt] = min(results);
        xopt=xmin+stepsize*(iopt-1);
        depth = depth - 1;
    end
    
end