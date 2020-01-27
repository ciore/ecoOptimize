classdef GCMMAFunctions
  methods(Static)
    
    %%
    function gcmma = init(xval,xnam,xmin,xmax,maxiter)
      % This sets up the GCMMA parameters and starting points
      [f0val,fval,df0dx,dfdx] = LEnOpFunctions.optFunctions(xval,xnam,true);
      m=size(fval,1);
      n=size(xval,1);
      epsimin = 0.0000001;
      xold1   = xval;
      xold2   = xval;
      low     = xmin;
      upp     = xmax;
      c       = 1000*ones(m,1);
      d       = 1*ones(m,1);
      a0      = 1;
      a       = 0*ones(m,1);
      raa0    = 0.01;
      raa     = 0.01*ones(m,1);
      raa0eps = 0.000001;
      raaeps  = 0.000001*ones(m,1);
      kkttol  = 0;
      innerit = 0;
      outeriter = 0;
      maxoutit  = maxiter;
      iout = [outeriter innerit];
      xout = [xval'];
      fout = [f0val];
      cout = [fval'];
      for name=who'
        if not(ismember(name,{'xval','xnam','xmin','xmax'}))
          eval(['gcmma.',name{:},'=',name{:},';'])
        end
      end
    end
    
    %%
    function [gcmma,xval] = run(gcmma,xval,xnam,xmin,xmax,plotlive)
      % This is a slightly modified version of GCMMA-MMA-code-1.5/gctoymain.m
      %
      addpath('../GCMMA-MMA-code-1.5')
      %
      names=fieldnames(gcmma);
      for i=1:numel(names)
        eval([names{i},'=gcmma.',names{i},';'])
      end
      %%%% begin optimisation
      %%%% If outeriter=0, the user should now calculate function values
      %%%% and gradients of the objective- and constraint functions at xval.
      %%%% The results should be put in f0val, df0dx, fval and dfdx:
      if outeriter < 0.5
        [f0val,fval,df0dx,dfdx] = LEnOpFunctions.optFunctions(xval,xnam,true);
        innerit=0;
        iout = [outeriter innerit];
        xout = [xval'];
        fout = [f0val];
        cout = [fval'];
      end
      %
      %%%% The outer iterations start:
      kktnorm = kkttol+10;
      outit = 0;
      while kktnorm > kkttol & outit < maxoutit
        outit   = outit+1;
        outeriter = outeriter+1;
        %%%% The parameters low, upp, raa0 and raa are calculated:
        [low,upp,raa0,raa] = ...
          asymp(outeriter,n,xval,xold1,xold2,xmin,xmax,low,upp, ...
          raa0,raa,raa0eps,raaeps,df0dx,dfdx);
        %%%% The GCMMA subproblem is solved at the point xval:
        [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp] = ...
          gcmmasub(m,n,outeriter,epsimin,xval,xmin,xmax,low,upp, ...
          raa0,raa,f0val,df0dx,fval,dfdx,a0,a,c,d);
        %%%% The user should now calculate function values (no gradients)
        %%%% of the objective- and constraint functions at the point xmma
        %%%% ( = the optimal solution of the subproblem).
        %%%% The results should be put in f0valnew and fvalnew.
        [f0valnew,fvalnew] = LEnOpFunctions.optFunctions(xmma,xnam,false);
        %%%% It is checked if the approximations are conservative:
        [conserv] = concheck(m,epsimin,f0app,f0valnew,fapp,fvalnew);
        %%%% While the approximations are non-conservative (conserv=0),
        %%%% repeated inner iterations are made:
        innerit=0;
        if conserv == 0
          while conserv == 0 & innerit <= 15
            innerit = innerit+1;
            %%%% New values on the parameters raa0 and raa are calculated:
            [raa0,raa] = ...
              raaupdate(xmma,xval,xmin,xmax,low,upp,f0valnew,fvalnew, ...
              f0app,fapp,raa0,raa,raa0eps,raaeps,epsimin);
            %%%% The GCMMA subproblem is solved with these new raa0 and raa:
            [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp] = ...
              gcmmasub(m,n,outeriter,epsimin,xval,xmin,xmax,low,upp, ...
              raa0,raa,f0val,df0dx,fval,dfdx,a0,a,c,d);
            %%%% The user should now calculate function values (no gradients)
            %%%% of the objective- and constraint functions at the point xmma
            %%%% ( = the optimal solution of the subproblem).
            %%%% The results should be put in f0valnew and fvalnew:
            [f0valnew,fvalnew] = LEnOpFunctions.optFunctions(xmma,xnam,false);
            %%%% It is checked if the approximations have become conservative:
            [conserv] = concheck(m,epsimin,f0app,f0valnew,fapp,fvalnew);
          end
        end
        %%%% No more inner iterations. Some vectors are updated:
        xold2 = xold1;
        xold1 = xval;
        xval  = xmma;
        %%%% The user should now calculate function values and gradients
        %%%% of the objective- and constraint functions at xval.
        %%%% The results should be put in f0val, df0dx, fval and dfdx:
        [f0val,fval,df0dx,dfdx] = LEnOpFunctions.optFunctions(xval,xnam,true);
        %%%% The residual vector of the KKT conditions is calculated:
        [residu,kktnorm,residumax] = ...
          kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s, ...
          xmin,xmax,df0dx,fval,dfdx,a0,a,c,d);
        iout = [iout; outeriter innerit];
        xout = [xout; xval'];
        fout = [fout; f0val];
        cout = [cout; fval'];
        for i=1:numel(names)
          eval(['gcmma.',names{i},'=',names{i},';'])
        end
        GCMMAFunctions.dispIter(gcmma)
        if plotlive
          GCMMAFunctions.plotIter(gcmma),drawnow
        end
      end
    end
    
    %%
    function plotIter(gcmma)
      subplot(3,2,1)
      plot(gcmma.iout(:,1),gcmma.fout)
      xlabel('iter'),ylabel('f0val')
      subplot(3,2,2)
      semilogy(gcmma.iout(2:end,1),abs(gcmma.fout(2:end,:)-gcmma.fout(1:end-1,:))./abs(gcmma.fout(1,:)))
      xlabel('iter'),ylabel('df0val/di')
      subplot(3,2,3)
      plot(gcmma.iout(:,1),gcmma.xout)
      xlabel('iter'),ylabel('xval')
      subplot(3,2,4)
      semilogy(gcmma.iout(2:end,1),abs(gcmma.xout(2:end,:)-gcmma.xout(1:end-1,:))./abs(gcmma.xout(1,:)))
      xlabel('iter'),ylabel('dxval/di')
      subplot(3,2,5)
      plot(gcmma.iout(:,1),gcmma.cout)
      xlabel('outit'),ylabel('fval')
      subplot(3,2,6)
      semilogy(gcmma.iout(2:end,1),abs(gcmma.cout(2:end,:)-gcmma.cout(1:end-1,:))./abs(gcmma.cout(1,:)))
      xlabel('iter'),ylabel('dfval/di')
    end
    
    %%
    function dispIter(gcmma)
      if gcmma.iout(end,1)==1
        fprintf(['GCMMA iterations:\n'])
        fprintf(['   iter   nini   f0val',repmat(' ',1,12-5),'xval',repmat(' ',1,gcmma.n*12-4),'fval-fmax\n'])
      end
      format=['%7i%7i',repmat('%12.3e',1,1+gcmma.n+gcmma.m),'\n'];
      fprintf(format,[gcmma.iout(end,:) gcmma.fout(end,:) gcmma.xout(end,:) gcmma.cout(end,:)])
    end
    
  end
end