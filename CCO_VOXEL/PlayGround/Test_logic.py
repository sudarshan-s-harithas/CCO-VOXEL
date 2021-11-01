import numpy as np
import mmd_cost as mmd
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import pandas as pd 
import csv 
from scipy import signal
num_samples_in_distribution =50
num_pts_sample_CEM = 50
iterations = 5
num_top_samples = 20
# fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})
impulse = signal.unit_impulse(35 , 17)
 #np.random.normal(0, 0.01, num_samples_in_distribution)
from mpl_axes_aligner import shift




# zero_dist =np.random.normal(0, 1 ,  num_samples_in_distribution)


noise = np.random.normal(0, 1 ,  num_samples_in_distribution)

def initilize():

	starting_point = np.array( [0,0])
	obstacle_location = np.array( [1,1])
	varx = 1
	vary = 1


	return starting_point , obstacle_location , varx, vary


starting_point , obstacle_location ,  init_varx , init_vary = initilize()


def perturb_point( mean_point , varx , vary ):

	perturb_setx  = np.random.normal(mean_point[0], varx, num_pts_sample_CEM)
	perturb_sety  = np.random.normal(mean_point[1], vary, num_pts_sample_CEM)

	perturb_set = np.zeros(( num_pts_sample_CEM , 2) )

	for i in range(num_pts_sample_CEM):
		

			perturb_set[i][0] = perturb_setx[i]
			perturb_set[i][1] = perturb_setx[i]

	return perturb_set



def draw_ellipse( point, covariance , iter_num):

	fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})


	if covariance.shape == (2, 2):
		U, s, Vt = np.linalg.svd(covariance)
		angle = np.degrees(np.arctan2(U[1, 0], U[0, 0]))
		width, height = 2 * np.sqrt(s)
	else:
		angle = 0
		width, height = 2 * np.sqrt(covariance)	


	nsig = 1

	width = covariance[0][0]
	height = covariance[1][1]

	v = np.arange( 0, 1, 0.3)


	for nsig in v:

		ellipse = Ellipse(point, nsig*width, nsig*height, angle=angle, alpha=0.25)
		ax.add_artist(ellipse)

	ax.set_xlim(-2.1 , 0)
	ax.set_ylim( 1, 0)

	plt.show()







# print( np.shape(perturb_set))

point = starting_point
varx = init_varx
vary = init_vary
ideal_distribution = np.zeros( (1, num_samples_in_distribution))
one_distribution = np.ones( (1, num_samples_in_distribution))
weights = np.ones( (1, num_samples_in_distribution) )



def get_distance_distribution(query_point):

	dist = np.sqrt( np.linalg.norm( query_point -obstacle_location )) 

	dist_distribution =  noise + dist  #zero_dist + dist
	dist_distribution = np.maximum( (0.75*one_distribution - dist_distribution), ideal_distribution) 

	return dist_distribution





def plot_pseudo_graph(distribution , iter_num):

	distribution += 0.1 
	# if( iter_num == 0):
	# 	distribution = new_dist 
	# mean  = np.round(  np.sum(distribution)/num_samples_in_distribution , decimals=3)
	# var = np.round(np.var(distribution) , decimals=3)

	impulse =  np.random.normal(0.6, 0.1, 50)

	
	s = pd.Series(distribution)*0.5

	#
	s2 = pd.Series(impulse)*0.5
	ax = s.plot.kde(linewidth = 2.2 , color ='navy')
	# ax = s2.plot.kde( label = "Target Distribution" )
	plt.xlim([-0.1, 2.5])

	shift.xaxis(ax, 0.1, 0.5, True)
	plt.xlabel("Constraint Violation (m)" , fontsize = 20)
	plt.ylabel("Probablity Density" , fontsize=20)
	plt.xticks(fontsize=10)


	# plt.text( 1 , 2.5 , 'Mean'+ str(mean) )
	# plt.text( 1 , 2.0 , 'Variance'+ str(var) )
	# plt.tight_layout()
	# plt.title('Distance Distribution')
	# plt.legend(fontsize=13)





	# plt.hist(distribution , bins= np.arange( -0.001, 2, 0.1) )
	plt.show()



for i in range(iterations):
	cost_list = []
	cost_arr = np.zeros( (1, num_samples_in_distribution))
	col =0 

	distance_dist = get_distance_distribution(point )

	plot_pseudo_graph( distance_dist[0] , 0)


	perturb_points = perturb_point( point , varx , vary )
	# print(np.shape(perturb_points))
	# print( perturb_points[5])
	cost_dictionary = {}
	elite_sample_Set = np.zeros( (num_top_samples+1,2)  )
	elite_sample_Setx = np.zeros( (num_top_samples+1)  )
	elite_sample_Sety = np.zeros( (num_top_samples+1)  )
	top_sample_counter  =0 

	for Query_points in perturb_points:
		

		distance_dist = get_distance_distribution(Query_points )

		# print( ideal_distribution[0])

		cost = mmd.mmd_cot_per_point(  distance_dist[0] , ideal_distribution  , weights , num_samples_in_distribution )

		cost_dict ={ str(Query_points) :cost }
		cost_dictionary.update(cost_dict)

		cost_arr[0][ col]  = cost[0][0]
		col +=1


	top_index = np.argpartition(cost_arr, num_top_samples)

	# print(top_index)

	for k in top_index[0] :
		# elite_sample_Set.append( perturb_points[k] )
		elite_sample_Set[top_sample_counter][0] = perturb_points[k][0]
		elite_sample_Set[top_sample_counter][1] = perturb_points[k][1]
		elite_sample_Setx[top_sample_counter] = perturb_points[k][0]
		elite_sample_Sety[top_sample_counter] = perturb_points[k][1]

		if( top_sample_counter >=num_top_samples ):

			break
		top_sample_counter +=1


	point = np.sum( elite_sample_Set , axis = 0) /top_sample_counter
	varx = np.var( elite_sample_Setx)
	vary = np.var( elite_sample_Sety)

	covar = np.array( [ [  varx , 0] , 
		[0 , vary]])

	# draw_ellipse( point, covar ,i )

	

	distance_dist = get_distance_distribution(point )

	cost = mmd.mmd_cot_per_point(  distance_dist[0] , ideal_distribution  , weights , num_samples_in_distribution )

	print( "Position of Current Update Point is " +  str(point)  )

	print( " MMD Cost at the given Point "  +  str(cost)) 

	
	print( "Updated Variance in X "  +  str(varx)) 
	print( "Updated Variance in Y " +  str(vary)) 

	print("----********---------")

	mean_distance_dist = get_distance_distribution(point)
	mean_distance_dist_transpose = np.zeros( (num_samples_in_distribution , 1))

	#plot_pseudo_graph( mean_distance_dist[0] , i)



	

	





		


		





























