 function [solution] = pendshooting(T, targetangle, N, solution)
    tic   
    global controls 
    
    controls.N     = N;           % Number of control nodes
    controls.T     = T;           % Duration of the motion
    controls.times = (0:N-1)'*controls.T/(N-1);
    controls.targetangle = targetangle;
    
    if (nargin == 4)    % use previous solution as an initial guess for control
        load('solution')
        oldresult = solution;
        oldtime   = solution.t;			
        controls.times = (0:controls.N-1)'*solution.t(end)/(controls.N-1);		% sample times for new optimization
        x0 = interp1(oldtime, oldresult.u, controls.times);  
    else
         x0 = zeros(N,1);     % initial guess for controls
    end
    
    options.cl = zeros(2,1);              % Lower bounds on constraints
    options.cu = zeros(2,1);              % Upper bounds on constraints
  
% Set the IPOPT options
    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.max_iter              = 1000;
    options.ipopt.linear_solver         = 'mumps';
    options.ipopt.constr_viol_tol       = 1e-5;
    options.ipopt.compl_inf_tol         = 1e-5;
    options.ipopt.mu_strategy           = 'adaptive';		
    options.ipopt.print_level           = 5;
    options.ipopt.bound_frac            = 0.001;
    options.ipopt.bound_push            = options.ipopt.bound_frac;
    options.ipopt.tol                   = 1e-4;

% The callback functions.
    funcs.objective                     = @objective;
    funcs.gradient                      = @gradient;
    funcs.constraints                   = @constraints;
    funcs.jacobian                      = @jacobian;
    funcs.jacobianstructure             = @jacobian;
  
% Run IPOPT.
  [x, info] = ipopt(x0,funcs,options);
  com.t = toc;

  % Run simulation and make plots
  y0 = [0;0];
  controls.x = x;
  [tt,yy]    = ode23(@odefunc, [0 controls.T], y0);
  uu         = interp1(controls.times, controls.x,tt); 
  
  figure(1)
  subplot(3,1,1)
  plot(tt,yy(:,1)*180/pi,'r','linewidth',2);
  ylabel('Position (degree)');
  
  subplot(3,1,2)
  plot(tt,yy(:,2)*180/pi,'b', 'linewidth', 2);
  ylabel('Velocity (deg/s)');
  
  subplot(3,1,3)
  plot(tt,uu, 'g', 'linewidth', 2);
  ylabel('Torque (N.m)');
  xlabel('Time (s)');
  	
  % make movie of the solution
	disp('Hit ENTER to generate animation...');
	pause
    L = 1;
    X = yy(:,1);
    avi = VideoWriter('pendSH.avi');
    avi.FrameRate = 12;
    open(avi);
	figure(2);
	clf;
	set(gcf,'Position',[5 100 650 650]);
	set(gcf, 'color', 'white');
	s = 1.5*L;
	for i=1:size(X)
		plot([-s s],[0 0],'k','LineWidth',2);
		hold on
		plot([0 L*cos(X(i)-pi/2)], [0 L*sin(X(i)-pi/2)],'b-o','LineWidth',2);
		axis('equal');
		axis('square');
		axis([-s s -s s]);
		title(['t = ' num2str(tt(i),'%8.3f')]);
		if (i==1)
			F = getframe(gca);
			frame = [1 1 size(F.cdata,2) size(F.cdata,1)];
		else
			F = getframe(gca,frame);
		end
		writeVideo(avi,F);
		drawnow;
		hold off;
	end
	close(avi);
% 	close(2);
    
% Check whether the problem solved or not!
  if info.status == 0;
  solution.status = 'Solved';
  else info.status ~=0;
  solution.status = 'Failed';
  end  
  
% save the solution file
  solution.x = yy;
  solution.t = tt;
  solution.u = uu;
  solution.n = N;
  solution.info = info;
  solution.comt = com.t;
  solution.costfunc    = sum(x.^2) * (controls.T/controls.N);
  solution.targetangle = solution.x(end,1)*180/pi;
  filename = 'solution.mat';
  save(filename,'solution')
  
% ------------------------------------------------------------------------
function f = objective (x)
   f = mean(x.^2);
% -------------------------------------------------------------------------
function g = gradient (x)
   global controls  
   g = 2*x/controls.N;
% ----------------------------------------------------------------------
function c = constraints(x)
   global  controls 
   y0 = [0;0];
   controls.x    = x;
   [~,yy]        = ode23(@odefunc, [0 controls.T],y0);
   c = yy(end,:) - [controls.targetangle,0];
   c = c';
%--------------------------------------------
function ydot = odefunc(t,y)
  global  controls
  d     = 1;
  m     = 1;
  g     = 9.81;
  u     = interp1(controls.times, controls.x, t);
  v     = y(2); 
  vdot  = -m.*g.*d.*sin(y(1))+ u;
  ydot  = [v ; vdot];
% ----------------------------------------------------------------------
function J = jacobian (x)
  global controls
  hh = 1e-2;
  x  = zeros(controls.N,1);
  c  = constraints(x);
  J  = zeros(numel(c), numel(x));
  for i = 1:numel(x);
      xsave  = x(i);
      x(i)   = x(i) + hh;
      chh    = constraints(x);
      J(:,i) = (chh-c)/hh;
      x(i)   = xsave;
  end
  J = sparse(J);