%% FUNCTION NAME: primalfep
%  Calculate $f_{\epsilon}(\rho)$ function, where $\epsilon$ value is
%  carefully chosen. 
%%

function [fval,realEpsilon] = primalfep(rho,keyMap,krausOperators,options)
    
    defaultoptions.epsilon = 0; % 0<=epsilon<=1/(e(d'-1)), epsilon is related to perturbed channel
    defaultoptions.perturbation = 1e-16; % a small value added to the minimum eigenvalue;
    if nargin == 4
        if ~isfield(options,'epsilon')
            options.epsilon = defaultoptions.epsilon;
        end
        if ~isfield(options,'perturbation')
            options.perturbation = defaultoptions.perturbation;
        end
    else 
        options = defaultoptions;
    end


    if nargin == 3 || isempty(krausOperators)
        % in the case that there is no post-selection (so no Kraus operator).
        dim = size(rho,1);
        eigMin = lambda_min(rho);
        realEpsilon = 0;
        
        
        if eigMin <= 0
            % if we need to do a perturbation.
            
            % first attempt is to use the epsilon value supplied by the user if reasonable.
            if options.epsilon > 0 && options.epsilon <= 1/(exp(1)*(dim-1))
                realEpsilon = options.epsilon;
                rhoPrime = (1-realEpsilon) * rho + realEpsilon * eye(dim)/dim;
                eigMin2 = lambda_min(rhoPrime);
                % check again the perturbed rho.
                if eigMin2 <=0
                    % now give up the user's choice and choose a reasonable
                    % epsilon value
                     realEpsilon = (-eigMin + options.perturbation) * dim;
                     % check again realEpsilon is in the range where the
                     % theorem 3 can apply
                     if realEpsilon < 0 || realEpsilon > 1/(exp(1)*(dim-1))
                        ME = MException('primalfep:badRho','Please check your rho from the first step calculation. It is not suitable for the second step calculation');
                     	throw(ME);
                     else
                         rhoPrime = (1-realEpsilon) * rho + realEpsilon * eye(dim)/dim;
                     end % end of if for validity of realEpsilon
                    
                end % end of check perturbed rho
                   
            else
                realEpsilon = (-eigMin + options.perturbation) * dim;
                % check realEpsilon is in the range where the
                % theorem 3 can apply
                if realEpsilon < 0 || realEpsilon > 1/(exp(1)*(dim-1))
                	ME = MException('primalfep:badRho','Please check your rho from the first step calculation. It is not suitable for the second step calculation');
                  	throw(ME);
                else
                    rhoPrime = (1-realEpsilon) * rho + realEpsilon * eye(dim)/dim;
                end
            end % end of doing perturbation but before checking again
            
            eigMin3 = lambda_min(rhoPrime);
           
            % check again the perturbed rho.
                
            if eigMin3<=0
            	ME = MException('primalfep:badRhobadEpsilon','Please check your rho and epsilon');
                throw(ME);
            end
            rho = rhoPrime;
        end % end of perturbation on rho
      
   
    
        zRho = 0;
        for jMapElm = 1 : numel(keyMap)
            zRho = zRho + keyMap{jMapElm} * rho * keyMap{jMapElm};
        end
        
        eigMin4 = lambda_min(zRho);
        % check if Z(rho) is valid
        if eigMin4 <= 0
            epsilon2 = (-eigMin4+options.perturbation) * dim;
            zRho = (1-epsilon2) * zRho + epsilon2 * eye(dim)/dim;
            eigMin5 = lambda_min(zRho);
           
            % check again the perturbed Z(rho).
                
            if eigMin5<=0
            	ME = MException('primalfep:badRhobadKeymap','Please check your rho and key map');
                throw(ME);
            end
            realEpsilon = max(realEpsilon, epsilon2);
        end
    
    
        fval = real(trace(rho * (logm(rho) - logm(zRho))));
    
    else
        % in the case that there is a post-selection map.
        gRho = krausFunc(rho,krausOperators);
        dimPrime = size(gRho,1);
        eigMin = lambda_min(gRho);
        realEpsilon =0;
        
        if eigMin <= 0
            % if we need to do a perturbation.
            
            % first attempt is to use the epsilon value supplied by the user if reasonable.
            if options.epsilon > 0 && options.epsilon <= 1/(exp(1)*(dimPrime-1))
                
                realEpsilon = options.epsilon;
                gRhoPrime = (1-realEpsilon) * gRho + realEpsilon * eye(dimPrime)/dimPrime;
                eigMin2 = lambda_min(gRhoPrime);
                % check again the perturbed rho.
                if eigMin2 <=0
                    % now give up the user's choice and choose a reasonable
                    % epsilon value
                     realEpsilon = (-eigMin + options.perturbation) * dimPrime;
                     % check again realEpsilon is in the range where the
                     % theorem 3 can apply
                     if realEpsilon < 0 || realEpsilon > 1/(exp(1)*(dimPrime-1))
                        ME = MException('primalfep:badRho','Please check your rho from the first step calculation. It is not suitable for the second step calculation');
                     	throw(ME);
                     else
                         gRhoPrime = (1-realEpsilon) * gRho + realEpsilon * eye(dimPrime)/dimPrime;
                     end
                    
                 end
                   
            else
                realEpsilon = (-eigMin + options.perturbation) * dimPrime;
                % check realEpsilon is in the range where the
                % theorem 3 can apply
                if realEpsilon < 0 || realEpsilon > 1/(exp(1)*(dimPrime-1))
                        ME = MException('primalfep:badRho','Please check your rho from the first step calculation. It is not suitable for the second step calculation');
                        throw(ME);
                else
                    gRhoPrime = (1-realEpsilon) * gRho + realEpsilon * eye(dimPrime)/dimPrime;
                end
            end
            
            eigMin3 = lambda_min(gRhoPrime);
           
            % check again the perturbed rho.
                
            if eigMin3<=0
            	ME = MException('primalfep:badRhobadEpsilon','Please check your rho and epsilon');
                throw(ME);
            end
            gRho = gRhoPrime;
        end
      
   
    
        zRho = 0;
        
        for jMapElm = 1 : numel(keyMap)
            zRho = zRho + keyMap{jMapElm} * gRho * keyMap{jMapElm};
        end
        
        eigMin4 = lambda_min(zRho);
        if eigMin4 <= 0
            
            epsilon2 = (-eigMin4+options.perturbation) * dimPrime;
            zRho = (1-epsilon2) * zRho + epsilon2 * eye(dimPrime)/dimPrime;
            eigMin5 = lambda_min(zRho);
           
            % check again the perturbed Z(G(rho)).
                
            if eigMin5<=0
            	ME = MException('primalfep:badRhobadKeymap','Please check your rho and key map');
                throw(ME);
            end
            
            realEpsilon = max(realEpsilon, epsilon2);
        end
    
    
        fval = real(trace(gRho * (logm(gRho) - logm(zRho))));
    end

end