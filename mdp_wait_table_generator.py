import numpy as np
import math
import random
import xlsxwriter
import mdptoolbox
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
import time


def num_state(num_bin):

#gives the num of states, the in put is a matrix that gives the num of for each state variable 
#in the order range, altitde , relative horizontal velocity, relative vertical velocity and 
#relative turn angle

    st=len(num_bin)
    prev_st=1
    for i in range(st):
        total_states=prev_st*num_bin[i];
        prev_st=total_states;
    return total_states+2
    


def vertex_update(x,y,z,x_rate,y_rate,A_rate,transition_time):

    X=x+x_rate*transition_time;
    Y=y+y_rate*transition_time;
    H=z+A_rate*transition_time;
    
    
    return (X,Y,H)
   

def overlap_cube(occ_cube,trgt_cube):

    x_length=max(min(trgt_cube[1],occ_cube[1])-max(trgt_cube[0],occ_cube[0]),0);
    y_length=max(min(trgt_cube[3],occ_cube[3])-max(trgt_cube[2],occ_cube[2]),0);
    z_length=max(min(trgt_cube[5],occ_cube[5])-max(trgt_cube[4],occ_cube[4]),0);
    
    new_volume=x_length*y_length*z_length   
    occ_volume=(occ_cube[1]-occ_cube[0])*(occ_cube[3]-occ_cube[2])*(occ_cube[5]-occ_cube[4]);
    over_lap_cube=(new_volume/occ_volume);

    if over_lap_cube<0:

        print(over_lap_cube)
        raise Exception('Invalid overlap!Check for valid points')

    
    return over_lap_cube
    
 
def setup_model (x_rng,y_rng,alt,rel_hv,rel_vv,turn_rate,num_bin):
    """Input: range of values for each state
    number of bin for each state
    this set up the bin and value of each bin, gives vertices of the distance"""

    #span of the total bin
    
    x_rng_span=x_rng[1]-x_rng[0]
    y_rng_span=y_rng[1]-y_rng[0]
    alt_span=alt[1]-alt[0]
    rel_hv_span=rel_hv[1]-rel_hv[0]
    rel_vv_span=rel_vv[1]-rel_vv[0]
    turn_rate_span=turn_rate[1]-turn_rate[0] #htis is turn rate dynamics from intruder
    
    #value of each bin
   
    bin_val_xrng=(x_rng_span/num_bin[0]) 
    bin_val_yrng=(y_rng_span/num_bin[1])
    bin_val_alt=alt_span/num_bin[2]
    bin_val_hv=(rel_hv_span/(num_bin[3]-1))
    bin_val_vv=(rel_vv_span/(num_bin[4]-1))
    bin_val_turn_rate=(turn_rate_span/(num_bin[5]-1))
    
    #create a matrix to fill in the value 
    x_rng_vertices=np.arange(x_rng[0],x_rng[1]+bin_val_xrng,bin_val_xrng)
    y_rng_vertices=np.arange(y_rng[0],y_rng[1]+bin_val_yrng,bin_val_yrng)
    alt_vertices=np.arange(alt[0],alt[1]+bin_val_alt,bin_val_alt)
    rel_hv_vertices=np.arange(rel_hv[0],rel_hv[1]+bin_val_hv,bin_val_hv)
    rel_vv_vertices=np.arange(rel_vv[0],rel_vv[1]+bin_val_vv,bin_val_vv)
    turn_rate_vertices=np.arange(turn_rate[0],turn_rate[1]+bin_val_turn_rate,bin_val_turn_rate)  #no change just tae the turn_dynamics of intruder as it is
    
    
    bin_val=abs(np.array([bin_val_xrng,bin_val_yrng,bin_val_alt,bin_val_hv,bin_val_vv,bin_val_turn_rate]))

    return (bin_val,x_rng_vertices,y_rng_vertices,alt_vertices,rel_hv_vertices,rel_vv_vertices,turn_rate_vertices)
    
    

def get_idx (state_num, num_bin):
    """
    Input: bin number of current state and number of bin of each state
    Output: index number of the state(integar)
    I assume for each X_dis, Y_dis changes first and all 4 state remain same. As all possible combination of x_dis and y_dis iterates altitude changes.
    Everytime alitude changes all the combination of x_dis and y_dis iterates over while other three states remain same. After all possible combination of 
    x_dis,y_dis and altitude change relative horizontal velocity changes and contitues like this for rest of the state"""

    #global indx
    
    #state_num includes six values x_rng,y_rng,height,h_v,turn_rate,v_v for current state bin number (i.e bin number can not be zero)
    #num_bin includes no of bins for six states in above sequence 
    
    #find multiplying factor for vertical velocity bin, find the idx of v_v and multiply each box dimension to items

    position_factor=state_num[0]-1  #this gives what has passed, like if position factor is 1, that means all the y_dis value is iterated for x_bin=1; and currently at x_bin=2
    
    position_state1=position_factor*num_bin[1]+state_num[1]
    position_factor2=state_num[2]-1;
    position_state=position_state1+num_bin[0]*num_bin[1]*position_factor2
    
    velocity_factor=state_num[3]-1   
    velocity_factor2=state_num[4]-1
    velocity_state=(num_bin[0]*num_bin[1]*num_bin[2]*num_bin[3]*velocity_factor)+(num_bin[0]*num_bin[1]*num_bin[2])*velocity_factor2
    
    basic=position_state+velocity_state
 
    turn_factor=state_num[5]-1
    turn_state=(num_bin[0]*num_bin[1]*num_bin[2]*num_bin[3]*num_bin[4])*turn_factor

    indx=np.int(basic+turn_state)-1 #since start with 0 index thus -1
    
    return(indx)

def get_state(s,num_bin):
    """
    Input:s is the index of state_matrix and num_bin is an array of bin value
          num_bin=[x,y,h,rel_hv,rel_vv,turn_rate]      
    Output:Numpy array of the six bin number for current state 
    reverse the encoding of getidx 
    """
    s_num=s+1    #this is the state number
    
    factor=s_num//(num_bin[0]*num_bin[1]) 
    remain_state=s_num%(num_bin[0]*num_bin[1])
  

    if remain_state==0:
    
        y_bin=num_bin[1]    
        x_bin=num_bin[0]
    else:  
    
        r_factor=remain_state//num_bin[1]
        r_remain=remain_state%num_bin[1]
        if r_remain==0:
            y_bin=num_bin[1]
            if r_factor==0:
                x_bin=remain_state
            else:
                x_bin=r_factor;
        else:
            x_bin=r_factor+1;
            y_bin=r_remain;
          
    t_factor=s_num//(num_bin[0]*num_bin[1]*num_bin[2]*num_bin[3]*num_bin[4])
    
    t_remain=s_num%(num_bin[0]*num_bin[1]*num_bin[2]*num_bin[3]*num_bin[4])
    
    
    if t_remain==0:
        t_bin=t_factor
    else:
        t_bin=t_factor+1;
        
    
    sh_factor=s_num//(num_bin[0]*num_bin[1]*num_bin[2])
    sh_remain=s_num%(num_bin[0]*num_bin[1]*num_bin[2])
    
    
    if sh_remain==0:
    
        sh_num=(num_bin[0]*num_bin[1]*num_bin[2])
    else:
        sh_num=sh_remain
  
    h_factor=sh_num//(num_bin[0]*num_bin[1])
    
    h_remain=sh_num%(num_bin[0]*num_bin[1])
    
    
    if h_remain==0:
    
        h_bin=h_factor
       
    else: 
    
        h_bin=h_factor+1
    
    
    if t_remain==0:
        sv_num=num_bin[0]*num_bin[1]*num_bin[2]*num_bin[3]*num_bin[4]
        
    else:
    
        sv_num=t_remain
        
    
    v_factor=sv_num//(num_bin[0]*num_bin[1]*num_bin[2]*num_bin[4])
    v_remain=sv_num%(num_bin[0]*num_bin[1]*num_bin[2]*num_bin[4])

    if v_remain==0:
    

        rv_bin=num_bin[4];
        rh_bin=v_factor 
       
        
    else:  
        
        rv_factor=v_remain//(num_bin[0]*num_bin[1]*num_bin[2])
        rv_remain=v_remain%(num_bin[0]*num_bin[1]*num_bin[2])

        if rv_remain==0:
            rv_bin=rv_factor
            
        else:

            rv_bin=rv_factor+1

        rh_bin=v_factor+1
    
    current_state=np.array([x_bin,y_bin,h_bin,rh_bin,rv_bin,t_bin])
    
    return current_state
    
  
    
def reward(hmd,current_dh,hmd_threshold,dh_value,time_factor,transition_time):
    """" generate reward matrix
    Input:hmd, current_dh, collision threshold values, factor, time interval
    Output: reward matrix"""

    R_h=hmd_threshold/hmd
    R_v=dh_value/current_dh
    R_d=-np.minimum(R_h,R_v)
    R_t=time_factor*transition_time
    R_w=R_d+R_t

    return R_w
    
def corrected_heading(heading_angle,turn_rate):
    """to correct the heading angle for intruder while taking a negative turn"""
    correct_heading=heading_angle+turn_rate;
    if correct_heading<0:
        correct_heading=360+correct_heading
  

    return correct_heading
    
    
def collision_parameter_gen(ownship_parameter, intruder_parameter,current_cube):
    """
    Verified with different values
    Input:
    Ownship parameters:  horizontal velocity, corrected_heading(including turn_rate)
    Intruder_parameters: horizontal velocity, corrected_heading(including turn_rate)
    relative_parameter:  current_cube
    Output:
    Current hmd value, tau mod and time to closest point of approach """


    intruder_velocity=intruder_parameter[0]
    intruder_corrected_heading=intruder_parameter[1]

    
    ownship_velocity=ownship_parameter[0]
    ownship_corrected_heading=ownship_parameter[1]

    
    xrng1=current_cube[0]
    xrng2=current_cube[1]
    yrng1=current_cube[2]
    yrng2=current_cube[3]
    zrng1=current_cube[4]
    zrng2=current_cube[5]
            
    if xrng1>=0:
        
        dx=0-xrng2   # dx=own_pos-int_pos ; for this model, ownship is always at 0; 
    else:

        dx=0-xrng1
    if yrng1>=0:
       
        dy=0-yrng2
    else:
 
        dy=0-yrng1
        
    rng=abs(math.sqrt(dx**2+dy**2))  
        
    v_rx=(ownship_velocity*math.cos(math.radians(ownship_corrected_heading)))-(intruder_velocity*math.cos(math.radians(intruder_corrected_heading)))
    v_ry=(ownship_velocity*math.sin(math.radians(ownship_corrected_heading)))-(intruder_velocity*math.sin(math.radians(intruder_corrected_heading)))

    
    R_rate=(dx*v_rx+dy*v_ry)/rng
  
    
    if rng==0:
        time_to_cpa=0;
    else:
        time_to_cpa=-(dx*v_rx+dy*v_ry)/(v_rx**2+v_ry**2)
        
    
    
    if time_to_cpa>=0:
        current_hmd=math.sqrt((dx+v_rx*time_to_cpa)**2+(dy+v_ry*time_to_cpa)**2)
        
    else:
        current_hmd=10000000  #just a big value 

    
    if rng<=dmod_threshold:
        tau_modified=0
        
    elif rng>dmod_threshold and R_rate<0:
        tau_modified=-((rng**2-dmod_threshold**2)/(rng*R_rate))
        
    else:
        tau_modified=10000000      #just a big value 
    
    if zrng2>0:
        current_dh_value=zrng2
    
    else:
        current_dh_value=abs(zrng1)
        
    return [current_hmd,current_dh_value,tau_modified,time_to_cpa]
    
def check_duplicate(new_cube,start_value,bin_val):
    """check if vertices go out of state and return unique sets of cube vertices where the points ends up"""
    x_points=[new_cube[0],new_cube[1]]
    y_points=[new_cube[2],new_cube[3]]
    z_points=[new_cube[4],new_cube[5]]
    x_start=start_value[0]
    y_start=start_value[1]
    z_start=start_value[2]
    
    set_cube=[]
    set_cube_vertices=[]
    points=[]
    out_bound_vertice=0
    collsion_bound=0
    for xx in range(len(x_points)):
        for yy in range(len(y_points)):
            for zz in range(len(z_points)):
                idx_x=int(abs(x_start-x_points[xx])//bin_val[0])  #these are indx
                idx_y=int(abs(y_start-y_points[yy])//bin_val[1])
                idx_z=int(abs(z_start-z_points[zz])//bin_val[2])
                
                if idx_x==num_bin[0]:
                
                    idx_x_state=idx_x  #these are bin number of state number
                else:
                    idx_x_state=idx_x+1
                    
                if idx_y==num_bin[1]:
                    idx_y_state=idx_y
                    
                else:
                    idx_y_state=idx_y+1
                    
                    
                if idx_z==num_bin[2]:
                    idx_z_state=idx_z
                else:
                    idx_z_state=idx_z+1
                    
                if abs(x_points[xx])>R or abs(y_points[yy])>R or abs(z_points[zz])>Z : 
                
                    out_bound_vertice=out_bound_vertice+1
     
                    continue
                    
                vertices_idx=[idx_x_state,idx_y_state,idx_z_state]
                point=[x_points[xx],y_points[yy],z_points[zz]]
                points.append(point)
                cube_vertices=[x_rng_vertices[idx_x_state-1],x_rng_vertices[idx_x_state],y_rng_vertices[idx_y_state-1],y_rng_vertices[idx_y_state],z_rng_vertices[idx_z_state-1],z_rng_vertices[idx_z_state]]
 
                set_cube.append(vertices_idx)
                set_cube_vertices.append(cube_vertices)
                
    unique_cube_idx_set = set(map(tuple,set_cube))  #need to convert the inner lists to tuples so they are hashable
    unique_cube_idx = list((map(list,unique_cube_idx_set)))

   
    return (unique_cube_idx,out_bound_vertice)
    
    
def matrix_filler(flattend_transition_matrix,state_num,target_idx,probability_value):
    """Insert the value of calculated probability in the target index of matrix"""
    flatten_idx=state_num*states+target_idx
 
    current_probability_value=flattend_transition_matrix[flatten_idx]
    updated_probability_value=current_probability_value+probability_value
    
    np.put(flattend_transition_matrix,flatten_idx,updated_probability_value) 
    
def boundary_dynamics(threshold_dynamics,state_dynamics):
    """We assumed aircraft can't exceed a turn_rate of 5 deg/s; so when intruder is already in a bin where it's taking a 5 deg/s turn 
    it cannot take a turn_rate more than 5 deg/s, this will prevent this occurance. Similarly, it ensures that no motion model cause a relative velocity
    to exceed the preset start and end value; THIS ASSUMPTION CAN BE CHANGED
    Input: Threshold values for relative horizontal velocity (min,max), relative vertical velocity (min,max) and turn_rate(negative,positive)
    current_dynamics current bin value for above sequence"""
    
    filter_dynamics=[];
 
    numrows=0;
    for row in threshold_dynamics:
        
        if state_dynamics[numrows]< row[0]: #less than min
            current_dynamics=row[0]
          
        elif state_dynamics[numrows]> row[1]:
            current_dynamics=row[1]
        else:
            current_dynamics=state_dynamics[numrows]
            
        filter_dynamics.append(current_dynamics)
        
        numrows+=1;

    return filter_dynamics
    


def head_angle_estimation(xy_ind):
    """Creates a 2-D array of angle for different vertices; assume x and y range is equal"""

    
    num_of_ind=len(xy_ind)
    max_value=max(xy_ind)

    heading_angle_matrix=np.zeros((num_of_ind,num_of_ind))
    heading_angle_matrix_flat=np.ndarray.flatten(heading_angle_matrix);
    num_of_rectange=int((max_value)/2)
    mid_point_x=(max_value//2)
    mid_point_y=(max_value//2)
    
    
    
    dict_square={};
    for xx in range(num_of_rectange):
        array_ind=np.arange(xx,max_value+1-xx)
        list_ind=list(array_ind)
        
        dict_square[xx]=list_ind

    
    for xx in range(num_of_rectange):
        mat_col=dict_square[xx]
        num_of_col=len(mat_col)
        min_ind=min(mat_col)
        max_ind=max(mat_col)

    
        for ii in range(num_of_col):
            indx_180=xx*num_of_ind+mat_col[ii]
            indx_0=(max_value-xx)*num_of_ind+mat_col[ii]
   
            
            if mat_col[ii]==min_ind:
                np.put(heading_angle_matrix_flat,indx_0,225)
                np.put(heading_angle_matrix_flat,indx_180,315)
            elif mat_col[ii]==max_ind:
                np.put(heading_angle_matrix_flat,indx_0,135)
                np.put(heading_angle_matrix_flat,indx_180,45)
            else:
                np.put(heading_angle_matrix_flat,indx_0,0)
                np.put(heading_angle_matrix_flat,indx_180,180)

            idx_to_fill=[i for i in xy_ind if (i>min_ind and i<max_ind)]
       
            for yy in range(len(idx_to_fill)):

                indx_side=idx_to_fill[yy]*num_of_ind+ii
                if ii>mid_point_x:
                    np.put(heading_angle_matrix_flat,indx_side,90)
                elif ii==mid_point_x:
                    np.put(heading_angle_matrix_flat,indx_side,0)
                else:
                    np.put(heading_angle_matrix_flat,indx_side,270)

        heading_angle_estimates=np.reshape(heading_angle_matrix_flat,(num_of_ind,num_of_ind))
    
    return heading_angle_estimates
    
def bin_to_idx(bin,max_num_bin):
    if bin==max_num_bin:
        indx=bin
    else:
        indx=bin+1
    return indx
    
def matrix_to_list(transition_matrix,whole_list=[]):

    """this funtion converts the transition matrix to a list of indices; for example if mdp is at state 1 it return the state where it ends up; like in state 4 5 6"""
    for xx in range(len(transition_matrix)):
        
        new_list=np.nonzero(transition_matrix[xx])
        new_list_ind=new_list[0];
        new_list_number=[n+1 for n in new_list_ind]
        whole_list.append(new_list_number)
        
  
    return whole_list
   
def list_to_adjacent(final_state_number,whole_list):

    """it return an adjaceny list ; for example, the state number (not the index of state, the integer state number) is the keys and 
    where it ends up is the list; this is just a dictionary form of the previous list it is needed for tracing all possible path"""
    
    list_dict={};
    for xx in range(final_state_number):
        state_key=xx+1;
        
        list_dict[state_key]=whole_list[xx]
    return list_dict
 


def locate_all_possible_paths(adjacency_list, start, end, path=[]):

    """this locates all the possible path starting and ending at certain node
    Input: Adjaceny list and start points are every state except from out state and collision state and end is collision state"""
    
    path = path + [start]
    if start == end:
        return [path]
    if start not in adjacency_list:
        return []
    paths = []
    for node in adjacency_list[start]:
        if node not in path:
            newpaths = locate_all_possible_paths(adjacency_list, node, end, path)
            for newpath in newpaths:
                paths.append(newpath)
    return paths
   
    
    

def calc_wait_time(all_paths,transition_matrix,transit_time):

    """this function calculate the wait time; the wait time is weighted with the probability value
    paths end at collision state; we first calculate the transition time of whole path and then is weighted with associted probability
    for example if the there's a transition from state 1 to collision state following 5 transitions, we calculated the probability of taking that path  by 
    multipying each state transition probability and then weighted the transition_time"""
    
    avg_wait_time=[];

    for row in all_paths:
        state_path=row

        num_of_path=len(state_path)
        if num_of_path==0:
        
            avg_wait_time.append(1000) #default value
            continue
          
        state_weighted_time=0;
        for row in state_path:
            num_of_state=len(row)
            ind=[n-1 for n in row]  #get the index to access the probability value form transition matrix_to_list
  

            total_probability=1;
            total_transition_time=0
            for xx in range(num_of_state-1): #no need to acces last one
                ii=ind[xx]
                jj=ind[xx+1]
                prob_value=transition_matrix[ii][jj];

                total_probability=total_probability*prob_value;
                total_transition_time+=transit_time
            
            weighted_total_path_time=total_probability*total_transition_time;
     
                
            state_weighted_time+=weighted_total_path_time          
    
        avg_state_wait_time=(state_weighted_time)  
        avg_wait_time.append(avg_state_wait_time)
   
    return avg_wait_time  
"""These are the inputs to code"""

"""Ownship Default value/Constant Value"""

own_pos=[0,0,0]   #ownship position
own_velocity_profile=np.array([30,5,0]) 

"""User_Input"""

##For the Intialization of model##
R=800 #this is range for both x and y
Z=200 
start_xrng=[-R,R] #smaller value first#
start_yrng=[-R,R]
start_alt=[-Z,Z]
start_hv=[80,200]
start_vv=[-10,10]
start_turn=[-5,5]
num_bin=[8,8,4,5,5,5]  #R_x,R_y,Z,Rhdot,Zdot and turnrate
transition_time=1; #1 second


##Intruder dynamic band##
#preset band; these values can be change#
intruder_horizontal_velocity_dynamics=np.array([[-10,-7.5,-5,-2.5,0,2.5,5,7.5,10,12],[0.05,0.05,0.10,0.25,0.10,0.15,0.10,0.10,0.05,0.05]]) #first row velocity value, second row probability#
intruder_vertical_velocity_dynamics=np.array([[-3,-2.5,-1,0,1,2.5,3],[0.10,0.15,0.15,0.20,0.15,0.15,0.10]])
intruder_turn_rate_dynamics=np.array([[-5,-2.5,0,2.5,5],[0.20,0.20,0.20,0.20,0.20]])
maneuver_probability=np.array([0.70,0.30])                                                       #first element horizontal,second element vertical#1
type_of_encounter=1 #head to head and overtaking if the value is 2 that is merging from left/right scenario

## Threshold values for collision in meter##
hmd_threshold=1220
dh_value_threshold=130
dmod_threshold=1220
tau_mod_threshold=35

"""Input ends"""

    
if __name__ == '__main__':
    
    #Intializing the model#    
    
    states=num_state(num_bin)
    horizontal_values=intruder_horizontal_velocity_dynamics[0][:]
    vertical_values=intruder_vertical_velocity_dynamics[0][:]
    turn_rate_values=intruder_turn_rate_dynamics[0][:]
    horizontal_values_ps=intruder_horizontal_velocity_dynamics[1][:]
    vertical_values_ps=intruder_vertical_velocity_dynamics[1][:]
    turn_rate_values_ps=intruder_turn_rate_dynamics[1][:]

    [bin_val,x_rng_vertices,y_rng_vertices,z_rng_vertices,rel_hv_bin_values,rel_vv_bin_values,turn_rate_bin_values]=setup_model(start_xrng,start_yrng,start_alt,start_hv,start_vv,start_turn,num_bin)
    print(len(x_rng_vertices))
    rng_indices=list(range(len(x_rng_vertices)))
    print(rng_indices)
    intduer_angle_vertices=head_angle_estimation(rng_indices)


    print(x_rng_vertices)
    print(y_rng_vertices)
    print(z_rng_vertices)
    print(rel_hv_bin_values)
    print(rel_vv_bin_values)
    print(turn_rate_bin_values)

    #Intialize the transfer matrix#

    transfer_mat_initial=np.zeros((states,states),dtype=np.float)
    transfer_mat_flat=np.ndarray.flatten(transfer_mat_initial)
    out_state=states-1
    col_state=states-2
    num_of_collision_state=0
    
    #Intialization of reward matrix#
    R_mat=np.zeros((states))
    
    for s in range(states-2):
        
        current_state=get_state(s,num_bin)
        x_bin=current_state[0]
        y_bin=current_state[1]
        z_bin=current_state[2]
        rel_hv_bin=current_state[3]
        rel_hv_idx=rel_hv_bin-1   #we need indx to access the value from array 
        rel_vv_bin=current_state[4]
        rel_vv_idx=rel_vv_bin-1
        turn_rate_bin=current_state[5]
        turn_rate_idx=turn_rate_bin-1
        
        xrng1=x_rng_vertices[x_bin-1]  
        xrng2=x_rng_vertices[x_bin]
        yrng1=y_rng_vertices[y_bin-1]
        yrng2=y_rng_vertices[y_bin]
        zrng1=z_rng_vertices[z_bin-1] 
        zrng2=z_rng_vertices[z_bin]
        current_cube=[xrng1,xrng2,yrng1,yrng2,zrng1,zrng2]
        
        intruder_horizontal_velocity=rel_hv_bin_values[rel_hv_idx]-own_velocity_profile[0]
        intruder_vertical_velocity=rel_vv_bin_values[rel_vv_idx]-own_velocity_profile[1]
        intruder_turn_rate=turn_rate_bin_values[turn_rate_idx]
        
        
        #intruder_heading=decide_intruder_heading(xrng2,yrng2,type=2)
        
        if xrng2<=0:
            col_ind_angle=x_bin-1
        else:
            col_ind_angle=x_bin
            
        if yrng2<=0:
            row_ind_angle=y_bin-1
        else:
            row_ind_angle=y_bin
            

        intruder_heading=intduer_angle_vertices[row_ind_angle,col_ind_angle]
        #corrected angle add here
        heading_corrected=corrected_heading(intruder_heading,intruder_turn_rate)
        
        
        intruder_param=[intruder_horizontal_velocity,intruder_vertical_velocity,heading_corrected]
        
        [current_hmd,current_dh_value,tau_mod,t_cpa]=collision_parameter_gen(own_velocity_profile,intruder_param,current_cube) 
        
        #Reward matrix#
        time_factor=1;
        R1=reward(current_hmd,current_dh_value,hmd_threshold,dh_value_threshold,time_factor,transition_time)
        np.put(R_mat,s,R1)

        if (current_hmd<=hmd_threshold and current_dh_value<=dh_value_threshold):# or (0<=tau_mod or tau_mod<=tau_mod_threshold):  #collision logic
            
            #if in collision states stays in collision states, no transition is calculated 
            
            num_of_collision_state+=1
            continue
        

        ##Transition to other states with horizontal dynamic model##
        horizontal_iter=0
        for horizontal_value in horizontal_values:
            turn_iter=0 #this i used to access the associated probability value from array
            for turn_rate_value in turn_rate_values:
            
                updated_state_dynamics=[rel_hv_bin_values[rel_hv_idx]+horizontal_value,rel_vv_bin_values[rel_vv_idx],turn_rate_bin_values[turn_rate_idx]+turn_rate_value]
                threshold_dynamics=[start_hv,start_vv,start_turn]
            
                filter_dynamics=boundary_dynamics(threshold_dynamics,updated_state_dynamics)
                
                updated_intruder_h_velocity=filter_dynamics[0]-own_velocity_profile[0]
                updated_intruder_v_velocity=filter_dynamics[1]-own_velocity_profile[1]
                updated_intruder_t_rate=filter_dynamics[2]-own_velocity_profile[2]
                
                heading_corrected=corrected_heading(intruder_heading,updated_intruder_t_rate)
                intruder_new_param=[updated_intruder_h_velocity,updated_intruder_v_velocity,updated_intruder_t_rate]
                
                intruder_hx=(updated_intruder_h_velocity)*math.cos(math.radians(heading_corrected))
                intruder_hy=(updated_intruder_h_velocity)*math.sin(math.radians(heading_corrected))
                own_hx=own_velocity_profile[0]*math.cos(math.radians(0))
                own_hy=own_velocity_profile[0]*math.sin(math.radians(0))
     
                rel_hx=intruder_hx-own_hx
                rel_hy=intruder_hy-own_hx
                z_rate=rel_vv_bin_values[rel_vv_idx]
                
                x_vertices=[xrng1,xrng2]
                y_vertices=[yrng1,yrng2]
                z_vertices=[zrng1,zrng2]
                
                
                new_vertices=[]
                new_x_vertices=[]
                new_y_vertices=[]
                new_z_vertices=[]
                
                for aa in range(len(x_vertices)):
                    for bb in range(len(y_vertices)):
                        for cc in range(len(z_vertices)):
                            new_vertice=vertex_update(x_vertices[aa],y_vertices[bb],z_vertices[cc],rel_hx,rel_hy,z_rate,transition_time)
                            new_x_vertice=new_vertice[0]
                            new_y_vertice=new_vertice[1]
                            new_z_vertice=new_vertice[2]
                            new_vertices.append(new_vertice)
                            new_x_vertices.append(new_x_vertice)
                            new_y_vertices.append(new_y_vertice)
                            new_z_vertices.append(new_z_vertice)
                            
                unique_x_vertices=sorted(set(new_x_vertices))
                unique_y_vertices=sorted(set(new_y_vertices))
                unique_z_vertices=sorted(set(new_z_vertices))
                
                new_occ_cube=unique_x_vertices
                new_occ_cube.extend(unique_y_vertices)
                new_occ_cube.extend(unique_z_vertices)
              
                
                out_cond=(unique_x_vertices[0]>=R or abs(unique_x_vertices[1])>=R or unique_y_vertices[0]>=R or abs(unique_y_vertices[1])>= R or unique_z_vertices[0]>=Z or abs(unique_z_vertices[1])>=Z)
                
                ##find the index of probability value##
                
                ii=horizontal_iter   
                jj=turn_iter
                
                if out_cond==True:   #whole matrix outside the range 
                
                    target_idx=out_state
                    p_value=(horizontal_values_ps[ii]*turn_rate_values_ps[jj]*maneuver_probability[0])

                    matrix_filler(transfer_mat_flat,s,target_idx,p_value)
                   
                    continue
                 
                #if not continue for rectangle 
                start_value=[start_hv[0],start_vv[0],start_turn[0]]
                [unique_cube_idx,out_bound_vertice]=check_duplicate(new_occ_cube,start_value,bin_val)
                
                new_vertices_set=[]
                for row in unique_cube_idx:
                    new_point=[x_rng_vertices[row[0]-1],x_rng_vertices[row[0]],y_rng_vertices[row[1]-1],y_rng_vertices[row[1]],z_rng_vertices[row[2]-1],z_rng_vertices[row[2]]]
                    new_vertices_set.append(new_point)
                    
                all_overlap=[]
                
                
                
                for row in new_vertices_set:
                    
                    overlap_value=overlap_cube(new_occ_cube,row)
                    all_overlap.append(overlap_value)
                
                total_overlap_value=round(sum(all_overlap),4)
               
                if total_overlap_value>1.000:  #less than 1 is okay because some points may end up outside 
                    print(total_overlap_value)
                    raise Exception('Somethings wrong with overlapping')

                #find the dynamics index#
                
                rel_hv_new_idx=(abs(start_hv[0]-filter_dynamics[0])//bin_val[3])
                rel_vv_new_idx=(abs(start_vv[0]-z_rate)//bin_val[4])
                turn_rate_new_idx=np.argwhere(turn_rate_bin_values==filter_dynamics[2])
                turn_rate_new_idx=turn_rate_new_idx[0]
                rel_hv_trgt_bin=bin_to_idx( rel_hv_new_idx, num_bin[3])
                rel_vv_trgt_bin=bin_to_idx(rel_vv_new_idx,num_bin[4])
                turn_rate_trgt_bin=turn_rate_new_idx+1
                
                dynamic_bin=[rel_hv_trgt_bin,rel_vv_trgt_bin,turn_rate_trgt_bin] #we need bin number of the state number while encoding the "intergaer state"
                
                
                p_value_filled=0 
                
                for aa in range(len(all_overlap)):
                

                    p_value=(horizontal_values_ps[ii]*turn_rate_values_ps[jj]*maneuver_probability[0]*all_overlap[aa])
                    pos_bin=unique_cube_idx[aa]
                    pos_bin.extend(dynamic_bin)
                    target_idx=get_idx(pos_bin,num_bin)
                 
                    p_value_filled+=p_value; #this indicates whether if every manuever ends up in some states 
                    matrix_filler(transfer_mat_flat,s,target_idx,p_value)
                    
                #now we find the value of probability that goes to out state (any probability that is not included in above mentioned states goes to out state
                
                final_p_value=p_value_filled
                row_p_value=horizontal_values_ps[ii]*turn_rate_values_ps[jj]*maneuver_probability[0]
                if final_p_value<row_p_value:
                    
                    out_fill=row_p_value-final_p_value
                    target_idx=out_state
             
                    matrix_filler(transfer_mat_flat,s,target_idx,out_fill)
                turn_iter+=1
            horizontal_iter+=1

        #for vertical maneuver#
        vertical_iter=0;
        for vertical_value in vertical_values:        
            
            
            updated_state_dynamics=[rel_hv_bin_values[rel_hv_idx],rel_vv_bin_values[rel_vv_idx]+vertical_value,turn_rate_bin_values[turn_rate_idx]] ;#adding the vertical change, everyother remain same
            threshold_dynamics=[start_hv,start_vv,start_turn]
        
            filter_dynamics=boundary_dynamics(threshold_dynamics,updated_state_dynamics)
            
            updated_intruder_h_velocity=filter_dynamics[0]-own_velocity_profile[0]
            updated_intruder_v_velocity=filter_dynamics[1]-own_velocity_profile[1]
            updated_intruder_t_rate=filter_dynamics[2]-own_velocity_profile[2]
            
            heading_corrected=corrected_heading(intruder_heading,updated_intruder_t_rate)
            
            intruder_hx=(updated_intruder_h_velocity)*math.cos(math.radians(heading_corrected))
            intruder_hy=(updated_intruder_h_velocity)*math.sin(math.radians(heading_corrected))
            own_hx=own_velocity_profile[0]*math.cos(math.radians(0))
            own_hy=own_velocity_profile[0]*math.sin(math.radians(0))
 
            rel_hx=intruder_hx-own_hx
            rel_hy=intruder_hy-own_hy
            z_rate=rel_vv_bin_values[rel_vv_idx]; #vertcial value added
            
            
            x_vertices=[xrng1,xrng2]
            y_vertices=[yrng1,yrng2]
            z_vertices=[zrng1,zrng2]
            
            
            new_vertices_v=[]
            new_x_vertices_v=[]
            new_y_vertices_v=[]
            new_z_vertices_v=[]
            
            for aa in range(len(x_vertices)):
                for bb in range(len(y_vertices)):
                    for cc in range(len(z_vertices)):
                        new_vertice_v=vertex_update(x_vertices[aa],y_vertices[bb],z_vertices[cc],rel_hx,rel_hy,z_rate,transition_time)
                        new_x_vertice_v=new_vertice_v[0]
                        new_y_vertice_v=new_vertice_v[1]
                        new_z_vertice_v=new_vertice_v[2]
                        new_vertices.append(new_vertice_v)
                        new_x_vertices_v.append(new_x_vertice_v)
                        new_y_vertices_v.append(new_y_vertice_v)
                        new_z_vertices_v.append(new_z_vertice_v)
                        
            unique_x_vertices_v=sorted(set(new_x_vertices_v))
            unique_y_vertices_v=sorted(set(new_y_vertices_v))
            unique_z_vertices_v=sorted(set(new_z_vertices_v))
            
            new_occ_cube_v=unique_x_vertices_v
            new_occ_cube_v.extend(unique_y_vertices_v)
            new_occ_cube_v.extend(unique_z_vertices_v)
          
            
            out_cond=(unique_x_vertices_v[0]>=R or abs(unique_x_vertices_v[1])>=R or unique_y_vertices_v[0]>=R or abs(unique_y_vertices_v[1])>= R or unique_z_vertices_v[0]>=Z or abs(unique_z_vertices_v[1])>=Z)
            
            ##find the index of probability value##
            
            kk=vertical_iter
            
            if out_cond==True:   #whole matrix outside the range 
            
                target_idx=out_state
                p_value_v=(vertical_values_ps[kk]*maneuver_probability[1])

                matrix_filler(transfer_mat_flat,s,target_idx,p_value_v)
               
                continue
             
            #if not continue for rectangle 
            start_value=[start_hv[0],start_vv[0],start_turn[0]]
            [unique_cube_idx_v,out_bound_vertice_v]=check_duplicate(new_occ_cube_v,start_value,bin_val);
            
            new_vertices_set_v=[]
            for row in unique_cube_idx_v:
                new_point=[x_rng_vertices[row[0]-1],x_rng_vertices[row[0]],y_rng_vertices[row[1]-1],y_rng_vertices[row[1]],z_rng_vertices[row[2]-1],z_rng_vertices[row[2]]]
                new_vertices_set_v.append(new_point)
                
            all_overlap_v=[]
            
            for row in new_vertices_set_v:
                
                overlap_value_v=overlap_cube(new_occ_cube,row)
                all_overlap_v.append(overlap_value_v)
            
            total_overlap_value_v=round(sum(all_overlap_v),4)
           
            if total_overlap_value_v>1.0:  

                raise Exception('Somethings wrong with overlapping')

            #find the dynamics index#
            
            rel_hv_new_idx=(abs(start_hv[0]-filter_dynamics[0])//bin_val[3])
            rel_vv_new_idx=(abs(start_vv[0]-filter_dynamics[1])//bin_val[4])
            turn_rate_new_idx=np.argwhere(turn_rate_bin_values==filter_dynamics[2])
            turn_rate_new_idx=turn_rate_new_idx[0]
            rel_hv_trgt_bin=bin_to_idx( rel_hv_new_idx, num_bin[3])
            rel_vv_trgt_bin=bin_to_idx(rel_vv_new_idx,num_bin[4])
            turn_rate_trgt_bin=turn_rate_new_idx+1
            
            dynamic_bin=[rel_hv_trgt_bin,rel_vv_trgt_bin,turn_rate_trgt_bin]  #we need bin number of the state number while encoding the "intergaer state"
            
            
            p_value_filled_v=0; 
            
            for aa in range(len(all_overlap)):
                p_value_v=(vertical_values_ps[kk]*maneuver_probability[1]*all_overlap[aa])
                pos_bin=unique_cube_idx[aa]   #bin number for position state 
                pos_bin.extend(dynamic_bin)
                target_idx=get_idx(pos_bin,num_bin)
        
                p_value_filled_v+=p_value_v; #this indicates whether if every manuever ends up in some states (including collision states)
                matrix_filler(transfer_mat_flat,s,target_idx,p_value_v)
                
            #now we find the value of probability that goes to out state (any probability that is not included in above mentioned states goes to out state
            final_p_value_v=p_value_filled_v  #sometimes it add up like 1.00001 mess up with a negative value letter to prevent that I rounded it up to 4 decimal places
            
            row_p_value_v=vertical_values_ps[kk]*maneuver_probability[1]
            if final_p_value_v<row_p_value_v: #if all the states end up in some other step doesn't go through this
                
                out_fill_v=row_p_value_v-final_p_value_v
                target_idx=out_state

                matrix_filler(transfer_mat_flat,s,target_idx,out_fill_v)
            
            vertical_iter+=1;
  
    trans_mat_shaped=np.reshape(transfer_mat_flat,(states,states)) 


    #reward matrix formatted as the toolbox's requred#
    C_mat=np.zeros((states)); #this the reward for taking evasive maneuevr action
    col_reward=-1000; #large negative reward at collision state
    np.put(R_mat,col_state,col_reward)
    np.put(C_mat,col_state,col_reward)
    
    
    states_where_less=np.argwhere(R_mat<0)
    states_to_delete=states_where_less[0:len(states_where_less)-1]; #don't want to delete the out state, where the reward is zero
    #np.savetxt('reward.csv',R_mat,delimiter=',')
    
    #np.savetxt('collision_state.csv',states_to_delete,delimiter=',')
    reward_action1_mat=np.delete(R_mat,states_to_delete)
    reward_action0_mat=np.delete(C_mat,states_to_delete) 
    
    reward_together=np.stack((reward_action0_mat,reward_action1_mat))
    
    final_reward_matrix=np.transpose(reward_together)

    #remove all the collision state and merged them to one#
    # I am deleting specific index with numpy delete since I had to find the index using np.where ; I get memory error if my array is too large
    #to avoid memory error I am batch processing 
    
    trans_mat_list=[]
    states_without_turn=num_bin[0]*num_bin[1]*num_bin[2]*num_bin[3]*num_bin[4]
    for i in range(len(turn_rate_bin_values)):
    
        batch_data=trans_mat_shaped[(i*states_without_turn):(i+1)*states_without_turn,:]
        row_to_del=[k for k in states_to_delete if (k<=(i+1)*states_without_turn and k>=(i*states_without_turn))] #rows that need to be deleted in this range
        batch_data_row=np.delete(batch_data,states_to_delete,axis=0)
        batch_data_col=np.delete(batch_data_row,states_to_delete,axis=1)
        trans_mat_list.extend(batch_data_col)
        
    #trans_mat_shaped_row=np.delete(trans_mat_shaped,states_to_delete,axis=0);
    #trans_mat_shaped_col=np.delete(trans_mat_shaped_row,states_to_delete,axis=1);
    list_col=[0]*len(trans_mat_list[0][:])
    trans_mat_list.append(list_col) 
    trans_mat_list.append(list_col)

    
    trans_mat_ready=np.array(trans_mat_list)

    final_state_number=len(trans_mat_ready[:][1]);
    
    new_out_state=final_state_number-1;
    new_col_state=final_state_number-2;
    row_probability_sum=np.sum(trans_mat_ready,axis=1);

    
    flat_trans_mat_final=np.ndarray.flatten(trans_mat_ready)
    
    
    for f in range(final_state_number-2):
        
        """if row_probability_sum[f]>1.0:
            print(f)
            print(row_probability_sum[f])
            raise Exception('sum is big somethings wrong')"""

        if row_probability_sum[f]<1.0: #there's collision state 
        
        
            collision_probability_value=1.0-row_probability_sum[f]
            flat_idx=f*final_state_number+new_col_state
            np.put(flat_trans_mat_final,flat_idx,collision_probability_value)

    flat_new_col_idx=new_col_state*final_state_number+new_col_state  #filling the outstate and collision as 1 i.e stays in coll/out states if in col/out state
    np.put(flat_trans_mat_final,flat_new_col_idx,1)
    flat_new_out_idx=new_out_state*final_state_number+new_out_state
    np.put(flat_trans_mat_final,flat_new_out_idx,1)
    
    #now prepare the transition matrix as per the toolbox requirement#
    """for the toolbox, we need to two transition matrix one for taking action 1 or wait another for taking 0 or evasive
    for taking evasive action the mdp will always ended up in out state if not in collision state"""
    
  
    transition_matrix_action1=np.reshape(flat_trans_mat_final,(final_state_number,final_state_number))
    print(np.sum(transition_matrix_action1))

    
    #np.savetxt('transition_matrix_action1.csv', transition_matrix_action1, delimiter=',') 
    
    #states of the wait#
    state_for_wait=(np.argwhere(R_mat>0)); #states that ar not collision or out state
    np.savetxt('wait_mat_states.csv',state_for_wait,delimiter=',')
    
    #states of the decision matrix#
    #for decision matrxi we will have decision in out and collision state too
    state_count_with_collision=np.concatenate((state_for_wait,np.array([states-1])),axis=None)
    state_for_decision=np.concatenate((state_count_with_collision,np.array([states])),axis=None)
    #np.savetxt('state_for_decision.csv',state_for_decision,delimiter=',')
    
    #the remaining code for solving the mdp and generating the decision matrix
    transfer_mat_for_maneuver=np.zeros((final_state_number,final_state_number),dtype=np.float)
    maneuver_mat=np.ndarray.flatten(transfer_mat_for_maneuver); #this is all zero matrix
    
    for x in range(final_state_number):
        if s==new_col_state:
            flat_idx=x*final_state_number+new_col_state
            np.put(maneuver_mat,flat_idx,1)
            
        else:
            flat_idx=x*final_state_number+new_out_state
            np.put(maneuver_mat,flat_idx,1)


    
    transition_matrix_action0=np.reshape(maneuver_mat,(final_state_number,final_state_number))
    
    transition_matrix_ready=np.stack((transition_matrix_action0,transition_matrix_action1))
    

    #solve with MDP
    vi = mdptoolbox.mdp.ValueIteration(transition_matrix_ready, final_reward_matrix,.99)
    vi.verbose
    vi.run()
    fun=vi.V;
    decision=np.array(vi.policy)
    

    #np.savetxt('decision_table.csv',decision,delimiter=',')
    print('Done MDP solving.. Starting to calculate wait time')
    
    #this creates the wait time#
    
    transition_matrix_remove_outstate=np.delete(transition_matrix_action1, -1, axis=1)
    transition_matrix_wait_calc=np.delete(transition_matrix_remove_outstate, -1, axis=0)
    total_state=len(transition_matrix_wait_calc);
    list_check=matrix_to_list(transition_matrix_wait_calc)
    get_list=list_to_adjacent(total_state,list_check)

    whole_path=[]
    for xx in range(total_state-1):

        #-1 because don't need collsion state#
        start_node=xx+1;
        
        end_node=total_state; #collision state

        path_for_state=locate_all_possible_paths(get_list,start_node,end_node)
     
        whole_path.append(path_for_state)

    wait_time=calc_wait_time(whole_path,transition_matrix_wait_calc,transit_time)
    
    np.savetxt('wait_time_table.csv',wait_time,delimiter=',')
    print(max(wait_time))