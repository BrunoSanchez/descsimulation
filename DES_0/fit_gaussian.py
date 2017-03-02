import numpy as np
from astropy.modeling import fitting,models
import random
import os
from astropy.nddata.utils import extract_array
from astropy.io import fits





def fit_gaussian(B,size):
	
	fitter=fitting.LevMarLSQFitter()
	y2,x2=np.mgrid[:size[0],:size[1]]
	
	p=models.Gaussian2D(x_mean=(size[0]/2.),y_mean=(size[1]/2.),x_stddev=1.,y_stddev=1.,theta=np.pi/4., 
	                    amplitude=B.max()-B.min())+models.Const2D(amplitude=B.min())
	
	out=fitter(p,x2,y2,B,maxiter=5000)
	return out


def extract_array_image(pix_array,x,y,fitshape,indices):

	position = (y, x)
	sub_array_data = extract_array(pix_array,
                                        fitshape,
                                        position,
                                        fill_value=-999.)
	return sub_array_data




	
def fit_gaussian_parche(iden,image,x,y,fitshape):
	
	pix_array=fits.getdata(image)
	
	indices=np.indices(pix_array.shape)
	
	param_out=np.zeros((len(x),8))
	
	for j in range(len(x)):
		estampa = extract_array_image(pix_array,x[j],y[j],fitshape,indices)
		out = fit_gaussian(estampa,fitshape)
		gauss_param=out[0]
		background=out[1]
		ymodel,xmodel=np.mgrid[:fitshape[0],:fitshape[1]]
		
		resid = np.sqrt(np.sum(np.square(estampa - out(xmodel,ymodel))))/(fitshape[0]*fitshape[1])
		
		elip = abs((gauss_param.x_stddev - gauss_param.y_stddev)/(gauss_param.x_stddev + gauss_param.y_stddev))
		th = np.rad2deg(gauss_param.theta)
		ab=gauss_param.x_stddev*gauss_param.y_stddev
		amp=gauss_param.amplitude.value
		if gauss_param.x_stddev > gauss_param.y_stddev:
			theta=(th- int(th/180)*180.)
		else:
			theta=((th- int(th/180)*180.)-90)
		if theta > 0:
			theta=180.-theta
		else:
			theta=abs(theta)
		if theta >180.:
			theta=theta-180.
		theta=180.-theta

		
		#~ import ipdb; ipdb.set_trace()		
		
		param_out[j,:]=np.array([iden[j],x[j],y[j],elip,theta,ab,amp,resid])
  
	return param_out        
