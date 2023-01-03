function a = InfillOpt()
%mdof[1,3]
% 1 total volume
% 3 upper bound


on = 1;
done = 0;
while on
    reset_data = strcat(pwd,filesep,'reset_data.mat');
    load(reset_data);
    if reset_flag == 1
        solve_flag = 0;
    end
    
    while solve_flag == 1
        % LOAD SETUP DATA
        % change this to your install directory
        setup_data = 'H:\www\data\InfillOpt\user_data.mat';
        load(setup_data);
        display('loaded setup data');

        % INITIALIZE X
        nelx = double(nelx);
        nely = double(nely);
        bVoid = bVoid.';
        bSolid = bSolid.';
        rmin = 1.6;
        nloop = 400;

        % Initialize internal parameters
        volfrac = 0.5;
        penal = 3;
        p = 8;
        r_hat = 6;
        move = 0.02;
        vol_min_inv_pNorm = (nelx*nely*(1-vol_min)^p)^(1/p);
        vol_max_pNorm = (nelx*nely*vol_max^p)^(1/p);

        %% MATERIAL PROPERTIES
        E0 = 1;
        Emin = 1e-9;
        nu = 0.3;

        %% PREPARE FINITE ELEMENT ANALYSIS
        A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
        A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
        B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
        B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
        KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
        nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
        edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
        edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
        iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
        jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);

        % DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
        Fsparse = sparse(2*(nely+1)*(nelx+1),1);
        Fsparse(Fx_ind,1) = Fx;
        Fsparse(Fy_ind,1) = Fy; 
        U = zeros(2*(nely+1)*(nelx+1),1);
        alldofs = [1:2*(nely+1)*(nelx+1)];
        freedofs = setdiff(alldofs,fixeddofs);


        % PREPARE FILTER
        iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
        jH = ones(size(iH));
        sH = zeros(size(iH));
        k = 0;

        for i1 = 1:nelx
          for j1 = 1:nely
              e1 = (i1-1)*nely+j1;
              if bVoid(j1,i1) == 0 && bSolid(j1,i1) == 0
                  for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
                      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                          if bVoid(j2,i2) == 0 && bSolid(j2,i2) == 0
                              e2 = (i2-1)*nely+j2;
                              k = k+1;
                              iH(k) = e1;
                              jH(k) = e2;
                              sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
                          end
                      end
                  end
              else
                  k = k+1;
                  iH(k) = e1;
                  jH(k) = e1;
                  sH(k) = rmin;
              end
          end
        end
        H = sparse(iH,jH,sH);
        Hs = sum(H,2);

        % INITIALIZE ITERATION
        x = repmat(volfrac,nely,nelx);
              for i1 = 1:nely
                  for j1 = 1:nelx
                      if bVoid(i1,j1) == 1
                          x(i1,j1) = 0;
                      end
                      if bSolid(i1,j1) == 1
                          x(i1,j1) = 1;
                      end
                  end
              end

        beta = 1;
        eta = 0.5;
        

        xTilde = x;
        xPhys = (tanh(beta*eta) + tanh(beta*(xTilde-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));

        loopbeta = 0;
        loop = 0;
        change = 1;

        xold1 = reshape(x,[nely*nelx,1]);
        xold2 = reshape(x,[nely*nelx,1]);
        low = 0;
        upp = 0;

        % START ITERATION
        while change > 0.001 && loop < nloop && solve_flag == 1
          loopbeta = loopbeta+1;
          loop = loop+1;
          % FE-ANALYSIS
          sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
          K = sparse(iK,jK,sK); K = (K+K')/2;
          U(freedofs) = K(freedofs,freedofs)\Fsparse(freedofs);

          % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
          ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
          c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
          dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;


          dv = ones(nely,nelx);
              [x_hat, x_hat_coef_neg, x_hat_coef] = average_hat_anisotropic(nelx,nely,r_hat,xPhys,p, bVoid);
              tmp = x_hat.^p;
              sum_x_hat_p = sum(sum(tmp));
              sum_x_hat_p_norm = sum_x_hat_p.^(1/p);

              x_inv_hat = 1-x_hat;
              tmp = x_inv_hat.^p;
              sum_x_inv_hat_p = sum(sum(tmp));
              sum_x_inv_hat_p_norm = sum_x_inv_hat_p.^(1/p);      

          % FILTERING/MODIFICATION OF SENSITIVITIES
            dx = beta * (1-tanh(beta*(xTilde-eta)).*tanh(beta*(xTilde-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));
            dc(:) = H*(dc(:).*dx(:)./Hs);
            dv(:) = H*(dv(:).*dx(:)./Hs);

          % UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES

        %     METHOD OF MOVING ASYMPTOTES (MMA)
              m = size(mdof,2);
              n = nelx*nely;  

              for i1 = 1:nely
                  for j1 = 1:nelx
                      if bVoid(i1,j1) == 1 || bSolid(i1,j1) == 1
                          dc(i1,j1) = 0;
                      end
                  end
              end
              df0dx = reshape(dc,[nelx*nely,1]);
              df0dx2 = 0*df0dx;
              dfdx = zeros(3,nelx*nely);
              dfdx(1,1:nelx*nely) = reshape(dv,[1,nelx*nely])/(nelx*nely*volfrac);
              dfdx(2,1:nelx*nely) = -sum_x_inv_hat_p(1,1,1).^(1/p-1) * reshape(x_hat_coef_neg(:,:,1), [nelx*nely,1]);
              dfdx(3,1:nelx*nely) =  sum_x_hat_p(1,1,1).^(1/p-1) * reshape(x_hat_coef(:,:,1), [nelx*nely,1]);
              dfdx(4,1:nelx*nely) = -sum_x_inv_hat_p(1,1,2).^(1/p-1) * reshape(x_hat_coef_neg(:,:,2), [nelx*nely,1]);
              dfdx(5,1:nelx*nely) = -sum_x_inv_hat_p(1,1,3).^(1/p-1) * reshape(x_hat_coef_neg(:,:,3), [nelx*nely,1]);
              dfdx(6,1:nelx*nely) =  sum_x_hat_p(1,1,2).^(1/p-1) * reshape(x_hat_coef(:,:,2), [nelx*nely,1]);
              dfdx(7,1:nelx*nely) =  sum_x_hat_p(1,1,3).^(1/p-1) * reshape(x_hat_coef(:,:,3), [nelx*nely,1]);

                  for ic = 2:7
                      tmp = reshape(dfdx(ic,:),[nely,nelx]);
                      for i1 = 1:nely
                          for j1 = 1:nelx
                              if bVoid(i1,j1) == 1 || bSolid(i1,j1) == 1
                                  tmp(i1,j1) = 0;
                              end
                          end
                      end
                      dfdx(ic,:) = reshape(H*(tmp(:).*dx(:)./Hs),[1,nelx*nely]);   
                  end

              dfdx2 = 0*dfdx;

              iter = loopbeta;
              xval = reshape(x,[nelx*nely,1]);
              xmin=max(0.0,xval-move);
              xmax=min(1,xval+move);

              f0val = c;
              fval = zeros(5,1);
              fval(1,1) = sum(sum(xPhys)) / (nelx*nely*volfrac) - 1;
              fval(2,1) = sum_x_inv_hat_p_norm(1,1,1) - vol_min_inv_pNorm;
              fval(3,1) = sum_x_hat_p_norm(1,1,1) - vol_max_pNorm;
              fval(4,1) = sum_x_inv_hat_p_norm(1,1,2) - vol_min_inv_pNorm;
              fval(5,1) = sum_x_inv_hat_p_norm(1,1,3) - vol_min_inv_pNorm;
              fval(6,1) = sum_x_hat_p_norm(1,1,2) - vol_max_pNorm;
              fval(7,1) = sum_x_hat_p_norm(1,1,3) - vol_max_pNorm;

              a0 = 1;
              a = zeros(m,1);
              c_ = ones(m,1)*1000;
              d = zeros(m,1);
              [xmma,ymma,zmma,lam,xsi,eta_,mu,zet,s,low,upp] = ...
              mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2,...
                  f0val,df0dx,df0dx2,fval(mdof),dfdx(mdof,:),dfdx2(mdof,:),low,upp,a0,a,c_,d);
              xnew = reshape(xmma,[nely,nelx]);
              xold2 = xold1;
              xold1 = xval;


              xTilde(:) = (H*xnew(:))./Hs;
              xPhys = (tanh(beta*eta) + tanh(beta*(xTilde-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));

              change = max(abs(xnew(:)-x(:)));
              x = xnew;


          % MODIFY BOUNDARY DENSITIES

              for i1 = 1:nely
                  for j1 = 1:nelx
                      if bVoid(i1,j1) == 1
                          xPhys(i1,j1) = 0;
                      end
                      if bSolid(i1,j1) == 1
                          xPhys(i1,j1) = 1;
                      end
                  end
              end

            % UPDATE RESULTS FILE
            xPhys_flipped = xPhys.';
            done = 0;
            SAVEPATH=strcat(pwd,filesep,'output');
            SAVEPATH2=strcat(pwd);

            if ( ~isdir(SAVEPATH))
                mkdir(SAVEPATH);
            end
            SAVEFILENAME=strcat(SAVEPATH,filesep,'result_data.mat');
            SAVEFILENAME2=strcat(SAVEPATH2,filesep,'reset_data.mat');

            disp(['Saved result to : ',SAVEFILENAME]);
            save(SAVEFILENAME, 'xPhys_flipped','done');

%             figure(1);
%             set(1, 'Position', [100, 450, 540, min(100+540*nely/nelx,540)]);
%               for i1 = 1:nely
%                   for j1 = 1:nelx
%                       if bVoid(i1,j1) == 1
%                           xPhys(i1,j1) = 0;
%                       end
%                       if bSolid(i1,j1) == 1
%                           xPhys(i1,j1) = 1;
%                       end
%                   end
%               end
%             colormap(gray); imagesc(-xPhys, [-1 0]); axis equal; axis tight; axis off; drawnow;

%             filename1 = sprintf('output\\rho-It%i.png',loop);
%             saveas(1,filename1,'png');
            
            % RESET
            load(SAVEFILENAME2);

            if reset_flag == 1
                reset_flag = 0;
                solve_flag = 0;
                xPhys_flipped = zeros(nelx,nely);
                done = 1;

                save(SAVEFILENAME, 'xPhys_flipped','done')
                save(SAVEFILENAME2,'solve_flag', 'reset_flag','stop_flag');
                display('Reset');

            elseif stop_flag == 1
                stop_flag = 0;
                solve_flag = 0;
                done = 1;
                save(SAVEFILENAME,'xPhys_flipped','done');
                save(SAVEFILENAME2,'reset_flag','solve_flag','stop_flag');
                display('Stopped');
       
            end

          % UPDATE HEAVISIDE REGULARIZATION PARAMETER
          if beta < 16 && (loopbeta >= 80 || change <= 0.002)
            beta = 2*beta;
            loopbeta = 0;
            change = 1;
            fprintf('Parameter beta increased to %g.\n',beta);
          end

        end

        % SAVE COMPLETION STATUS
        done = 1;
        save(SAVEFILENAME,'xPhys_flipped','done');
        display('Done');


    end
end



   

%%%%%%%%% AVERAGE FILTER - ANISOTROPIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_hat, x_hat_coef_neg, x_hat_coef]=average_hat_anisotropic(nelx,nely,rmin,x,p,bVoid)
x_hat=zeros(nely,nelx,3);
x_hat_coef_neg=zeros(nely,nelx,3);
x_hat_coef=zeros(nely,nelx,3);
for i = 1:nelx
  for j = 1:nely
    sum=0.0; 

        for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
          for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
              if (k-i)*(k-i) + (l-j)*(l-j) <= rmin*rmin
                      sum = sum+1;
                      x_hat(j,i,1) = x_hat(j,i,1) + x(l,k);
              end
          end
        end

    x_hat(j,i,1) = x_hat(j,i,1)/sum;
    

       for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
          for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
              if (k-i)*(k-i) + (l-j)*(l-j) <= rmin*rmin
                  if bVoid(l,k) == 0
                      x_hat_coef_neg(l,k,1) = x_hat_coef_neg(l,k,1) + (1-x_hat(j,i,1))^(p-1) * 1/sum;
                      x_hat_coef(l,k,1) = x_hat_coef(l,k,1) + (x_hat(j,i,1))^(p-1) * 1/sum;
                  end
              end
          end
       end
  end
end

for i = 1:nelx
  for j = 1:nely
    sum=0.0; 
    for k = max(i-floor(rmin/3),1):min(i+floor(rmin/3),nelx)
      for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
          if 9*(k-i)*(k-i) + (l-j)*(l-j) <= rmin*rmin
              sum = sum+1;
              x_hat(j,i,2) = x_hat(j,i,2) + x(l,k);
          end
      end
    end
    x_hat(j,i,2) = x_hat(j,i,2)/sum;
    
    for k = max(i-floor(rmin/3),1):min(i+floor(rmin/3),nelx)
      for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
          if 9*(k-i)*(k-i) + (l-j)*(l-j) <= rmin*rmin
              x_hat_coef_neg(l,k,2) = x_hat_coef_neg(l,k,2) + (1-x_hat(j,i,2))^(p-1) * 1/sum;
              x_hat_coef(l,k,2) = x_hat_coef(l,k,2) + (x_hat(j,i,2))^(p-1) * 1/sum;
          end
      end
    end
  end
end

for i = 1:nelx
  for j = 1:nely
    sum=0.0; 
    for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
      for l = max(j-floor(rmin/3),1):min(j+floor(rmin/3),nely)
          if (k-i)*(k-i) + 9*(l-j)*(l-j) <= rmin*rmin
              sum = sum+1;
              x_hat(j,i,3) = x_hat(j,i,3) + x(l,k);
          end
      end
    end
    x_hat(j,i,3) = x_hat(j,i,3)/sum;
    
    for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
      for l = max(j-floor(rmin/3),1):min(j+floor(rmin/3),nely)
          if (k-i)*(k-i) + 9*(l-j)*(l-j) <= rmin*rmin
              x_hat_coef_neg(l,k,3) = x_hat_coef_neg(l,k,3) + (1-x_hat(j,i,3))^(p-1) * 1/sum;
              x_hat_coef(l,k,3) = x_hat_coef(l,k,3) + (x_hat(j,i,3))^(p-1) * 1/sum;
          end
      end
    end
  end
end
