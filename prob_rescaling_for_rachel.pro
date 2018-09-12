

function predfunc_fiducial,x

	mu=2.44+0.00			
	sigma=0.28+0.00			
	
	y=(x-mu)/sigma
	y=0.5*(1.0-erf(y/sqrt(2.0)))
	return,y
	
end

pro rescale_probabilities_for_rachel,epsc,delta,cat,newp

	s=alog(cat.chisq/cat.chisq_mod)		
	bad=where(finite(s) eq 0,nbad)
	if nbad gt 0 then s[bad]=0.0

	makeaxis,x,-6.0,8.0,1000
	dof=4.0
	chisq=exp(x)
	y=1.0/2.0^(dof/2.0)/gamma(dof/2.0)*chisq^(dof/2.0)*exp(-chisq/2.0)  ;; chisq distribution
		
	xvec=alog(cat.chisq_mod)
	bad=where(cat.chisq_mod le exp(-6.0),nbad)
	if nbad gt 0 then xvec[bad]=-6.0
	xvec1=xvec+s
				
	rho=exp(-s)*interpol(y,x,xvec)
	rho0=interpol(y,x,xvec1)
		
	bad=where(rho lt 0.0,nbad)
	if nbad gt 0 then rho[bad]=min(y)
	bad=where(rho0 lt 0.0,nbad)
	if nbad gt 0 then rho0[bad]=min(y)
		
	gal_epschi=rho/rho0-1.0

	bad=where(gal_epschi ge 1e4,nbad)
	if nbad gt 0 then gal_epschi[bad]=1e4
		
	fblue=1.0-predfunc_fiducial(alog(cat.chisq_mod))
	gal_epsblue=fblue/(1.0-fblue)
	bad=where(1.0-fblue lt 1e-5,nbad)
	if nbad gt 0 then gal_epsblue[bad]=1e5	
			
	b=delta+gal_epschi+delta*gal_epschi
	a=gal_epsblue+epsc

	newp=cat.p*(1.0+b)/(1.0+cat.p*(a+b+a*b))
	bad=where(cat.r eq 0,nbad)
	if nbad gt 0 then newp[bad]=1.0				
			
end



pro remap_s82_chisq_to_dr8_chisq,chisq

	common remap,x,xeq
	
	if n_elements(x) eq 0 then begin	

		numbins=long(1e3)	
		makeaxis,x,0.0,30.0,numbins	

		dof=4.0
		y=x^(0.5*dof-1.0)*exp(-0.5*x)/gamma(0.5*dof)/2.0^(0.5*dof)	
		dof=3.0
		y3=x^(0.5*dof-1.0)*exp(-0.5*x)/gamma(0.5*dof)/2.0^(0.5*dof)
		cdf4=0.0*y
		cdf3=cdf4
		for i=1L,numbins-1 do begin
			cdf4[i]=int_tabulated(x[0:i],y[0:i])
			cdf3[i]=int_tabulated(x[0:i],y3[0:i])
		endfor
		xeq=interpol(x,cdf4,cdf3)
	endif
	
	chisq=interpol(xeq,x,chisq)

end



pro rescale_probabilities_for_rachel_s82,epsc,delta,cat,newp

	s=alog(cat.chisq/cat.chisq_mod)		
	bad=where(finite(s) eq 0,nbad)
	if nbad gt 0 then s[bad]=0.0

	makeaxis,x,-6.0,8.0,1000
	dof=3.0									;; use dof=3 for s82 because there is no u-band
	chisq=exp(x)
	y=1.0/2.0^(dof/2.0)/gamma(dof/2.0)*chisq^(dof/2.0)*exp(-chisq/2.0)  ;; chisq distribution
		
	xvec=alog(cat.chisq_mod)
	bad=where(cat.chisq_mod le exp(-6.0),nbad)
	if nbad gt 0 then xvec[bad]=-6.0
	xvec1=xvec+s
				
	rho=exp(-s)*interpol(y,x,xvec)
	rho0=interpol(y,x,xvec1)
		
	bad=where(rho lt 0.0,nbad)
	if nbad gt 0 then rho[bad]=min(y)
	bad=where(rho0 lt 0.0,nbad)
	if nbad gt 0 then rho0[bad]=min(y)
		
	gal_epschi=rho/rho0-1.0

	bad=where(gal_epschi ge 1e4,nbad)
	if nbad gt 0 then gal_epschi[bad]=1e4
	
	chisq=cat.chisq_mod
	remap_s82_chisq_to_dr8_chisq,chisq			;; turn a dof=3 chisq into a dof=4 chisq
	
	fblue=1.0-predfunc_fiducial(alog(chisq))
	gal_epsblue=fblue/(1.0-fblue)
	bad=where(1.0-fblue lt 1e-5,nbad)
	if nbad gt 0 then gal_epsblue[bad]=1e5
			
	b=delta+gal_epschi+delta*gal_epschi
	a=gal_epsblue+epsc

	newp=cat.p*(1.0+b)/(1.0+cat.p*(a+b+a*b))
				
			
end




pro test

	epsc=0.062
	delta=-0.095

	rmdir='/Users/erozo/Documents/Research/Clusters/redMaPPer/Catalogs/dr8_v5.10/'
	rmfile='dr8_run_redmapper_v5.10_lgt20_catalog_members_mod.fit'
	cat=mrdfits(rmdir+rmfile,1)
			
	zmin=0.1
	zmax=0.3	
	use=where(cat.z ge zmin and cat.z le zmax)
	cat=cat[use]
	
	rescale_probabilities,epsc,delta,cat,newp1
	rescale_probabilities_for_rachel,epsc,delta,cat,newp2
	
	plot,newp1,newp2,psym=3

end

