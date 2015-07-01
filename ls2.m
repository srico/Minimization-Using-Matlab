diary ('output.txt');
xc = fscanf(fopen('x.txt','r'),'%f');%Initial x input vector, to be read from a file
iter = 0 ;
fxc = f2(xc) ;%call to f2(), F vector function
minfc = minif(fxc) ;% call to minif(), f scalar function
minfn = minfc ;
a = 0.0001 ;%value of alpha
bt = 0 ;
ctr = 1 ;
lam0 = 1 ;%value of lambda
fprintf('\niter   xc(1)               xc(2)             F(1)               F(2)                f             lambda\n');		
fprintf('----------------------------------------------------------------------------------------------------------\n');
fprintf('%1d  % 10.8e   % 10.8e  % 10.8e   % 10.8e   % 10.8e   % 8.6f\n',ctr, xc(1), xc(2), fxc(1), fxc(2),minfc,lam0);
while (minfc > 1e-10)%main while loop; while f function value 
                     % is not in range, continue iterations
	iter = iter + 1 ;
	jxc = jaco (xc) ;% call to jaco(), jacobian function
	au = [jxc,-fxc] ;%creating augmented matrix
	rr = rref(au) ;%solving augmented matrix
	sn = [ rr(1,3) rr(2,3) ]' ;% Newton step
	xn = xc + sn ;% new value of x
	fxn = f2(xn) ;% new value of F
	minfn = minif(fxn);% new value of f
	lam0 = 1 ;% setting current value of lambda to 1
	while (minfn >= minfc -a*lam0*2*(minfc))% backtracking while loop
		bt = bt + 1 ; % counter of backtracking iterations
		lam1 = minfc/(minfn-minfc+2*minfc);%calculate lambda
		if (lam1 < 0.1*lam0 ) % choosing closer lambda boundary
			lam1 = 0.1*lam0 ; 
		elseif (lam1 > 0.5*lam0 )
			lam1 = 0.5*lam0;
		end
		xn = xc + lam1*sn ; 
		fxn = f2(xn);
		minfn = minif(fxn); 
		lam0 = lam1; 
	end% end of backtracking while loop
	xc = xn ; 
	fxc = fxn ;
	minfc = minfn ;
	ctr = ctr + 1 ;	% update counter	
fprintf('%1d  % 10.8e   % 10.8e  % 10.8e   % 10.8e   % 10.8e   % 8.6f\n',ctr, xc(1), xc(2), fxc(1), fxc(2),minfc,lam0);
end% end of main while loop
diary off ;







