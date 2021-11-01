import numpy as np
from numpy.linalg import multi_dot
import time 
from geometry_msgs.msg import PoseStamped
import rospy 
Query_point_publisher = rospy.Publisher('/QueryPoint', PoseStamped, queue_size=1)
from matplotlib import pyplot as plt 

EDT_distance =0 


def get_reduced_set(  result , desired_dist  , weights  , num_samples_of_distance_distribution ):

	extract_lower_bound = 45
	extract_upper_bound =55
	reduced_size = extract_upper_bound - extract_lower_bound

	reduced_actual_dist = result[extract_lower_bound:extract_upper_bound] 
	reduced_desired_dist = desired_dist[extract_lower_bound:extract_upper_bound] 


	kzx = np.zeros( ( reduced_size ,num_samples_of_distance_distribution )) 
	kz = np.zeros((reduced_size , reduced_size  ))

	for i in range(reduced_size):
		for j in range(reduced_size):

			kz[i][j] = (  reduced_actual_dist[i]*reduced_actual_dist[j] +1 )*(  reduced_actual_dist[i]*reduced_actual_dist[j] +1 )

	for i in range(reduced_size):
		for j in range(num_samples_of_distance_distribution):
			kzx[i][j] = (  reduced_actual_dist[i]*result[j] +1 )*(  reduced_actual_dist[i]*result[j] +1 )


	temp1 = (np.linalg.pinv(kz))
	##print( np.shape(temp1) , np.shape(kzx))
	temp2 = np.matmul( temp1 ,kzx )

	weights = np.matmul( temp2 , weights.T)
	weights =weights/100

	print("**************")

	print( np.sum(weights))
	print("**************")

	return  reduced_actual_dist , reduced_desired_dist, weights.T, reduced_size






def mmd_cot_per_point(  result , desired_dist  , weights , num_samples_of_distance_distribution ):

	# result, desired_dist , weights , num_samples_of_distance_distribution = get_reduced_set(  result , desired_dist  , weights , num_samples_of_distance_distribution)

	

	kernel_matrix1  = np.zeros( (num_samples_of_distance_distribution ,num_samples_of_distance_distribution) )


	for i in range(num_samples_of_distance_distribution):
		for j in range(num_samples_of_distance_distribution):

			kernel_matrix1[i][j] = (  result[i]*result[j] +1 )*(  result[i]*result[j] +1 )

	res1  = multi_dot([weights, kernel_matrix1  , np.transpose(  weights)])
	'''

	kernel_matrix2  = np.zeros( (num_samples_of_distance_distribution,num_samples_of_distance_distribution) )

	for i in range(num_samples_of_distance_distribution):
		for j in range(num_samples_of_distance_distribution):
					kernel_matrix2[i][j] =  (  result[i]*desired_dist[j] +1 )*(  result[i]*desired_dist[j] +1 )


	res2 = multi_dot([weights, kernel_matrix2  , np.transpose(  weights)])

	kernel_matrix3  = np.zeros( (num_samples_of_distance_distribution,num_samples_of_distance_distribution) )

	for i in range(num_samples_of_distance_distribution):
		for j in range(num_samples_of_distance_distribution):
			kernel_matrix3[i][j] = (  desired_dist[i]*desired_dist[j] +1 )*(  desired_dist[i]*desired_dist[j] +1 )

	res3  = multi_dot([weights, kernel_matrix3  , np.transpose(  weights)])
	'''

	res = res1 - num_samples_of_distance_distribution*num_samples_of_distance_distribution #2*res2 +res3

	return res



def EDT_callback(dist):

    global EDT_distance 

    EDT_distance = dist.data 
    # print(EDT_distance)
    is_edt_updated = 1





def mmd_cost_per_trajectory(trajectory_points , points_sampled , f):

	num_samples_of_distance_distribution = 100

	trajectory_cost = 0
	r_safe= 0.75

	global EDT_distance

	union_of_samples= [] 
	union_of_samples = np.array(union_of_samples)


	for i in range(points_sampled):
		x_point = trajectory_points[i][0]
		y_point = trajectory_points[i][1]
		z_point = trajectory_points[i][2]

		robot_position = np.array([ x_point , y_point , z_point] )
		##diff =obstacle_location - robot_position

		##d = np.sqrt( np.sum(   diff*diff))

		Query_point= PoseStamped()

		Query_point.pose.position.x = x_point 
		Query_point.pose.position.y = y_point 
		Query_point.pose.position.z = z_point 

		
		Query_point_publisher.publish(Query_point) 
		
		##print(Query_point)
		time.sleep(0.7)

		
		d = EDT_distance

		thresold_val = r_safe - ( d -2)  

		if( thresold_val > 0):
			# print(d, "entered the MMD "  + str(d) )

			d_distribution =  np.random.normal(d, 1.5,num_samples_of_distance_distribution)
			actual_distribution =  r_safe - d_distribution
			actual_distribution = np.maximum(  np.zeros( np.shape(actual_distribution) ), actual_distribution)
			union_of_samples =  np.append( union_of_samples  , actual_distribution)
			desired_dist = np.zeros( np.shape(actual_distribution))
			##ck.plot_probablity_densities( actual_distribution , desired_dist) 
			weights = np.ones( (1,num_samples_of_distance_distribution) )
			point_cost = mmd_cot_per_point( actual_distribution , desired_dist  , weights , num_samples_of_distance_distribution)
			trajectory_cost += point_cost


	f.write( str(union_of_samples))
	f.close()



	return trajectory_cost




def get_cost_per_trajectory(data_matrix , points_sampled):


	dim = np.shape(data_matrix)
	cost = np.zeros(( 1, dim[0]))



	for i in range( dim[0]):
		name = str(i) +"_traj.txt"

		f = open( name , "a")

		cost[0][i] = mmd_cost_per_trajectory( data_matrix[i][:][:] , points_sampled , f)
		print( "trajectory " + str(i) + " complete")

	min_index = np.argmin(cost)	

	return cost , min_index