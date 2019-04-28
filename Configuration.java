/* Copyright 2013, Gurobi Optimization, Inc. */

/* This example formulates and solves the following cloud assignment model:

	min	\sum_{j \in P} z_j a_j 
	s.t.	x_{ij} <= z_j, i \in V, j \in P
		y_{ikjl} <= x_{ij}, i \in V, j \in P, k \in R_i, l \in D_j
		\sum_{j \in P} \sum_{l \in D_j} y_{ikjl} = 1, i \in V, k \in R_i
		\sum_{j \in P} x_{ij} = 1, i \in V
		\sum_{k \in R_i} y_{ikjl} <= 1, i \in V, j \in P, l \in D_j
		\sum_{i \in V} \sum_{k \in R_i} v_{ik} y_{ikjl} <= S_{jl}, j \in P, l \in D_j
		\sum_{i \in V} \alpha_i x_{ij} <= C_j, j \in P
		\sum_{i \in V} \beta_i x_{ij} <= M_j, j \in P
		x_{ij}, y_{ikjl}, z_j binary
*/

import java.io.*;
import gurobi.*;

public class Configuration
{
	// data members
	int num_VM_type;		// number of VM types
	String [] VM_type;  		// VM types
	int [] VM_CPU;			// # VM CPUs
	double [] VM_memory; 		// # VM memory (GB)
	int [] VM_storage_num;		// # VM storage disk
	int [] VM_storage_size;		// # VM storage disk size (GB)

	int num_PM_type;		// number of PM types
	String [] PM_type;  		// PM types
	int [] PM_CPU;			// # VM CPUs
	double [] PM_memory;	 	// # VM memory (GB)
	int [] PM_storage_num;		// # VM storage disk
	int [] PM_storage_size;		// # VM storage disk size (GB)
	int [] PM_purchase_cost;	// # fixed cost (normalized)
	int [] PM_operation_cost;	// # fixed cost (normalized)

	int max_K;			// max num of VM storage disks for any VM Type
	int max_L;			// max num of PM storage disks for any PM Type

	//configuration
	int [] configuration;		// # a vector, each item is the number of VMs of type i assigned to a PM

	//the possible max assignment of each VM-PM pair
	int [] VM_PM_bounds;

	// constructors
	public Configuration(String conf_VM, String conf_PM)
	{
		String line = null;
		max_K = 0;
		max_L = 0;
 
		// wrap a BufferedReader around FileReader
		BufferedReader bufferedReader = null; 
 
		//Read VM configuration
		try{
			bufferedReader = new BufferedReader(new FileReader(conf_VM));

			//#VM_types
			line = bufferedReader.readLine();
			num_VM_type = Integer.parseInt(line);

			VM_type = new String[num_VM_type + 1];
			VM_CPU = new int[num_VM_type + 1];
			VM_memory = new double[num_VM_type + 1];
			VM_storage_num = new int[num_VM_type + 1];
			VM_storage_size = new int[num_VM_type + 1];

			//VM Type
			line = bufferedReader.readLine();

			for(int i = 1; i <= num_VM_type; i++)
			{
				line = bufferedReader.readLine();
				VM_type[i] = line;
			}

			//vCPU
			line = bufferedReader.readLine();

			for(int i = 1; i <= num_VM_type; i++)
			{
				line = bufferedReader.readLine();
				VM_CPU[i] = Integer.parseInt(line);
			}

			//memory (GiB)
			line = bufferedReader.readLine();

			for(int i = 1; i <= num_VM_type; i++)
			{
				line = bufferedReader.readLine();
				VM_memory[i] = Double.parseDouble(line);
			}

			//storage (GiB)
			line = bufferedReader.readLine();

			for(int i = 1; i <= num_VM_type; i++)
			{
				line = bufferedReader.readLine();
				int pos = line.indexOf("x");
				VM_storage_num[i] = Integer.parseInt(line.substring(0, pos - 1));
				if(max_K < VM_storage_num[i])
					max_K = VM_storage_num[i];
				VM_storage_size[i] = Integer.parseInt(line.substring(pos + 2));
			}
			bufferedReader.close();

		}catch(IOException e){
			System.err.println("Error: " + e.toString());
		}

		//Read PM configuration
		try{
			bufferedReader = new BufferedReader(new FileReader(conf_PM));

			//#PM_types
			line = bufferedReader.readLine();
			num_PM_type = Integer.parseInt(line);

			PM_type = new String[num_PM_type + 1];
			PM_CPU = new int[num_PM_type + 1];
			PM_memory = new double[num_PM_type + 1];
			PM_storage_num = new int[num_PM_type + 1];
			PM_storage_size = new int[num_PM_type + 1];
			PM_purchase_cost = new int[num_PM_type + 1];
			PM_operation_cost = new int[num_PM_type + 1];
		
			//PM Type
			line = bufferedReader.readLine();

			for(int i = 1; i <= num_PM_type; i++)
			{
				line = bufferedReader.readLine();
				PM_type[i] = line;
			}

			//vCPU
			line = bufferedReader.readLine();

			for(int i = 1; i <= num_PM_type; i++)
			{
				line = bufferedReader.readLine();
				PM_CPU[i] = Integer.parseInt(line);
			}

			//memory (GiB)
			line = bufferedReader.readLine();

			for(int i = 1; i <= num_PM_type; i++)
			{
				line = bufferedReader.readLine();
				PM_memory[i] = Double.parseDouble(line);
			}

			//storage (GiB)
			line = bufferedReader.readLine();

			for(int i = 1; i <= num_PM_type; i++)
			{
				line = bufferedReader.readLine();
				int pos = line.indexOf("x");
				PM_storage_num[i] = Integer.parseInt(line.substring(0, pos - 1));
				if(max_L < PM_storage_num[i])
					max_L = PM_storage_num[i];
				PM_storage_size[i] = Integer.parseInt(line.substring(pos + 2));
			}

			//purchase cost (normalized)
			line = bufferedReader.readLine();

			for(int i = 1; i <= num_PM_type; i++)
			{
				line = bufferedReader.readLine();
				PM_purchase_cost[i] = Integer.parseInt(line);
			}

			//operation cost (normalized)
			line = bufferedReader.readLine();

			for(int i = 1; i <= num_PM_type; i++)
			{
				line = bufferedReader.readLine();
				PM_operation_cost[i] = Integer.parseInt(line);
			}

			bufferedReader.close();

		}catch(IOException e){
			System.err.println("Error: " + e.toString());
		}
		
		configuration = new int[num_VM_type + 1];
		VM_PM_bounds = new int[num_VM_type + 1];
	}

	// methods

	// input: VM index; output: VM type index
	private int getVMTypeIndex(int VM_index)
	{
		int VM_Type_index = 0;

		int VM_total = 0;		
		for(int i = 1; i <= num_VM_type; i++)
		{
			VM_total += configuration[i];
			if(VM_total >= VM_index)
			{
				VM_Type_index = i;
				break;
			}
		}

		if(0 == VM_Type_index)
		{
			System.err.println("Error: wrong VM_index!");
			System.exit(1);
		}

		return VM_Type_index;
	}

	// check the validity of a VM-PM assignment 
	private boolean check_validity(int PM_type_index)
	{
		if(PM_type_index < 1 || PM_type_index > num_PM_type)
		{
			System.err.println("Error: wrong PM_type_index = " + PM_type_index);
			return false;
		}

		//check CPU assignment
		int N = 0;
		int CPU_amount = 0;
		for(int i = 1; i <= num_VM_type; i++)
		{
			N += configuration[i];
			CPU_amount += VM_CPU[i] * configuration[i];
		}

		if(CPU_amount > PM_CPU[PM_type_index])
		{
			System.err.println("Invalid configuration, too many CPU allocated.");
			return false;
		}

		//check memory assignment

		int memory_amount = 0;
		for(int i = 1; i <= num_VM_type; i++)
		{
			memory_amount += VM_memory[i] * configuration[i];
		}

		if(memory_amount > PM_memory[PM_type_index])
		{
			System.err.println("Invalid configuration, too many memory allocated.");
			return false;
		}
	
		//check disk assignment
		if(check_disk_assignment(N, PM_type_index) == false)
		{
			System.err.println("Invalid configuration, infeasible disk assigment.");
			return false;
		}

		return true;
	}

	// only need to check the feasibility of disk assignment, we solve an optimization problem instead
	/*
	min	\sum_{i \in V} \sum_{k \in R_i} \sum_{l \in D_j} \sum y_{ikl}  
	s.t.	\sum_{l \in D_j} y_{ikjl} = 1, i \in V, k \in R_i
		\sum_{k \in R_i} y_{ikl} <= 1, i \in V, l \in D_j
		\sum_{i \in V} \sum_{k \in R_i} v_{ik} y_{ikl} <= S_{l}, l \in D_j
		y_{ikl} binary
	*/
	private boolean check_disk_assignment(int N, int PM_type_index)
	{
		try {
			GRBEnv	 env	= new GRBEnv("disk_assignment.log");
			GRBModel  model = new GRBModel(env);

			// Create variables

			// y_{ikl}: y[i][k][l], i = 1 to N, k 
			GRBVar [][][] y = new GRBVar[N + 1][max_K + 1][max_L + 1];
			for(int i = 1; i <= N; i++)			
			{
				int VM_type_index = getVMTypeIndex(i);
				for(int k = 1; k <= VM_storage_num[VM_type_index]; k++)
				{
					for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
						y[i][k][l] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "y_{" + i + "," + k + "," + l + "}");
				}
			}
			
			// Integrate new variables
			model.update();

			// Set objective: min	\sum_{i \in V} \sum_{k \in R_i} \sum_{l \in D_j} \sum y_{ikl}  
			GRBLinExpr expr = new GRBLinExpr();
			for(int i = 1; i <= N; i++)			
			{
				int VM_type_index = getVMTypeIndex(i);
				for(int k = 1; k <= VM_storage_num[VM_type_index]; k++)
				{
					for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
						expr.addTerm(1.0, y[i][k][l]);
				}
			}
			model.setObjective(expr, GRB.MINIMIZE);

			// Add constraint c_c: 	\sum_{l \in D_j} y_{ikl} = 1, i \in V, k \in R_i
			for(int i = 1; i <= N; i++)			
			{
				int VM_type_index = getVMTypeIndex(i);
				for(int k = 1; k <= VM_storage_num[VM_type_index]; k++)
				{
					expr = new GRBLinExpr();
					for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
						expr.addTerm(1.0, y[i][k][l]);
					model.addConstr(expr, GRB.EQUAL, 1.0, "c_c_{" + i + "," + k + "}");
				}
			}

			// Add constraint c_e:	\sum_{k \in R_i} y_{ikl} <= 1, i \in V, l \in D_j
			for(int i = 1; i <= N; i++)			
			{
				for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
				{	
					expr = new GRBLinExpr();
					int VM_type_index = getVMTypeIndex(i);
					for(int k = 1; k <= VM_storage_num[VM_type_index]; k++)
						expr.addTerm(1.0, y[i][k][l]);
					model.addConstr(expr, GRB.LESS_EQUAL, 1.0, "c_e_{" + i + ","  + l + "}");
				}
			}

			// Add constraint c_f:	\sum_{i \in V} \sum_{k \in R_i} v_{ik} y_{ikl} <= S_{l}, j \in P, l \in D_j
			for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
			{
				expr = new GRBLinExpr();
				for(int i = 1; i <= N; i++)			
				{				
					int VM_type_index = getVMTypeIndex(i);
					for(int k = 1; k <= VM_storage_num[VM_type_index]; k++)
						expr.addTerm(VM_storage_size[VM_type_index], y[i][k][l]);
				}
				model.addConstr(expr, GRB.LESS_EQUAL, PM_storage_size[PM_type_index], "c_f_{" + l + "}");
			}

			// Update model

			model.update();

			// write the model to a file
			//model.write(output_file);

			// Optimize model

			model.optimize();

			// Print objective
			System.out.println("Obj: " + model.get(GRB.DoubleAttr.ObjVal));

			// Print solution
			for(int i = 1; i <= N; i++)			
			{
				int VM_type_index = getVMTypeIndex(i);
				for(int k = 1; k <= VM_storage_num[VM_type_index]; k++)
				{
					for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
						System.out.println(y[i][k][l].get(GRB.StringAttr.VarName)
								 + " " + y[i][k][l].get(GRB.DoubleAttr.X));
				}
			}

			// Dispose of model and environment

			model.dispose();
			env.dispose();

			return true;
		 } catch (GRBException e) {
			System.err.println("Error code: " + e.getErrorCode() + ". " + e.getMessage());
			return false;
		 }
	}

	// pre_compute the max number of assignments for each VM-PM pair
	public void pre_compute(int PM_type_index)
	{
		System.out.print("PM_type = " + PM_type_index + ": ");
		for(int i = 1; i <= num_VM_type; i++)
		{
			if(VM_storage_num[i] > PM_storage_num[PM_type_index] || VM_storage_size[i] > PM_storage_size[PM_type_index])
			{				
				VM_PM_bounds[i] = 0;

				if(i < num_VM_type)
					System.out.print(VM_PM_bounds[i] + ", ");
				else
					System.out.println(VM_PM_bounds[i]);

				continue;
			}

			// CPU
			int num = (int) (PM_CPU[PM_type_index] / VM_CPU[i]);
			VM_PM_bounds[i] = num;

			//if(VM_PM_bounds[i] > 8)
				//VM_PM_bounds[i] = 0;
			
			// memory
			if(PM_memory[PM_type_index] / VM_memory[i] < VM_PM_bounds[i])
				VM_PM_bounds[i] = (int) (PM_memory[PM_type_index] / VM_memory[i]);

			//one physical disk can have at most one virtual disk for one VM
			num = (int) (PM_storage_num[PM_type_index] / VM_storage_num[i]);
			if(num > 0)	
				num = (int)(PM_storage_size[PM_type_index] / VM_storage_size[i]);
		
			if(num < VM_PM_bounds[i])
				VM_PM_bounds[i] = num;

			if(i < num_VM_type)
				System.out.print(VM_PM_bounds[i] + ", ");
			else
				System.out.println(VM_PM_bounds[i]);
		}
	}


	private void enumerate_configurations(int VM_type_index, int PM_type_index, BufferedWriter output_writer)
	{
		for(int num = VM_PM_bounds[VM_type_index]; num >= 0; num--)
		{	
			configuration[VM_type_index] = num;
			if(VM_type_index == num_VM_type)
			{
				boolean bvalid = check_validity(PM_type_index);
				if(bvalid)
				{
					//it is a valid configuration, write it to the list file
					try{
						for(int i = 1; i <= num_VM_type - 1; i++)
							output_writer.write(configuration[i] + ",");
						output_writer.write(configuration[num_VM_type] + "\n");
				        	output_writer.flush();    
					}catch(IOException e){
						System.err.println("Error: " + e.toString());
					}

					//return;
				}/*else{
					for(int i = 1; i <= num_VM_type - 1; i++)
						System.out.print(configuration[i] + ",");
					System.out.println(configuration[num_VM_type] + "\n");
				}*/
			}else{
				enumerate_configurations(VM_type_index + 1, PM_type_index, output_writer);
			}
		}				
	}

	public void enumerate_configurations(int PM_type_index, String conf_output)
	{
		try{
			BufferedWriter output_writer = new BufferedWriter(new FileWriter(conf_output));
			enumerate_configurations(1, PM_type_index, output_writer);
			output_writer.close();
		}catch(IOException e){
			System.err.println("Error: " + e.toString());
		}
	}

	private void compute_resource_assignment(String conf_file, String output_file, int PM_type_index)
	{
		String line = null;

		//Read PM-VM configuration file
		try{
			BufferedReader bufferedReader = new BufferedReader(new FileReader(conf_file));
			BufferedWriter output_writer = new BufferedWriter(new FileWriter(output_file));

			//#configurations
			line = bufferedReader.readLine();
			int num_confs = Integer.parseInt(line);

			for(int j = 1; j <= num_confs; j++)
			{
				line = bufferedReader.readLine();
				int num_VMs = 0;
				int num_CPUs = 0;
				double used_memory = 0;
				int pre_pos = 0;
				int pos = 0;
				for(int i = 1; i <= num_VM_type - 1; i++)
				{
					pos = line.indexOf(",", pre_pos);
					configuration[i] = Integer.parseInt(line.substring(pre_pos, pos));
					pre_pos = pos + 1;
					num_VMs += configuration[i];
					num_CPUs += configuration[i] * VM_CPU[i];
					used_memory += configuration[i] * VM_memory[i];
				}
				configuration[num_VM_type] = Integer.parseInt(line.substring(pre_pos));
				num_VMs += configuration[num_VM_type];
				num_CPUs += configuration[num_VM_type] * VM_CPU[num_VM_type];
				used_memory += configuration[num_VM_type] * VM_memory[num_VM_type];

				//print "#VMs,#CPUs,used memory,used disk num,used disk size"				
				String resourceline = "";
				resourceline = resourceline + num_VMs + "," + num_CPUs + "," + used_memory + "," +	
compute_disk_assignment(num_VMs, PM_type_index);

				output_writer.write(resourceline + "\n");
			}
			bufferedReader.close();
			output_writer.close();

		}catch(IOException e){
			System.err.println("Error: " + e.toString());
		}

	}

	/*
	min	\sum_{l \in D} z_l  
	s.t.	y_{ikl} <= z_l, l \in D
		\sum_{l \in D} y_{ikl} = 1, i \in V, k \in R_i
		\sum_{k \in R_i} y_{ikl} <= 1, i \in V, l \in D
		\sum_{i \in V} \sum_{k \in R_i} v_{ik} y_{ikl} <= S_{l}, l \in D
		y_{ikl} z_l binary
	*/
	private String compute_disk_assignment(int N, int PM_type_index)
	{
		try {
			GRBEnv	 env	= new GRBEnv("disk_assignment.log");
			GRBModel  model = new GRBModel(env);

			// Create variables

			// y_{ikl}: y[i][k][l], i = 1 to N, k 
			GRBVar [][][] y = new GRBVar[N + 1][max_K + 1][max_L + 1];
			for(int i = 1; i <= N; i++)			
			{
				int VM_type_index = getVMTypeIndex(i);
				for(int k = 1; k <= VM_storage_num[VM_type_index]; k++)
				{
					for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
						y[i][k][l] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "y_{" + i + "," + k + "," + l + "}");
				}
			}

			// z_l: z[l], l \in D
			GRBVar [] z = new GRBVar[PM_storage_num[PM_type_index] + 1];
			for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)			
				z[l] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "z_{" + l + "}");
			
			// Integrate new variables
			model.update();

			// Set objective: min	\sum_{l \in D_j} z_{l}  
			GRBLinExpr expr = new GRBLinExpr();
			for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
			{
				expr.addTerm(1.0, z[l]);
			}
			model.setObjective(expr, GRB.MINIMIZE);

			// Add constraint c_a: 	y_{ikl} \leq z_l, l \in D
			for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
			{
				for(int i = 1; i <= N; i++)			
				{
					int VM_type_index = getVMTypeIndex(i);
					for(int k = 1; k <= VM_storage_num[VM_type_index]; k++)
					{
						expr = new GRBLinExpr();
						expr.addTerm(1.0, y[i][k][l]);
						expr.addTerm(-1.0, z[l]);
						model.addConstr(expr, GRB.LESS_EQUAL, 0.0, "c_a_{" + i + "," + k + "," + l + "}");
					}
				}
			}

			// Add constraint c_c: 	\sum_{l \in D_j} y_{ikl} = 1, i \in V, k \in R_i
			for(int i = 1; i <= N; i++)			
			{
				int VM_type_index = getVMTypeIndex(i);
				for(int k = 1; k <= VM_storage_num[VM_type_index]; k++)
				{
					expr = new GRBLinExpr();
					for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
						expr.addTerm(1.0, y[i][k][l]);
					model.addConstr(expr, GRB.EQUAL, 1.0, "c_c_{" + i + "," + k + "}");
				}
			}

			// Add constraint c_e:	\sum_{k \in R_i} y_{ikl} <= 1, i \in V, l \in D_j
			for(int i = 1; i <= N; i++)			
			{
				for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
				{	
					expr = new GRBLinExpr();
					int VM_type_index = getVMTypeIndex(i);
					for(int k = 1; k <= VM_storage_num[VM_type_index]; k++)
						expr.addTerm(1.0, y[i][k][l]);
					model.addConstr(expr, GRB.LESS_EQUAL, 1.0, "c_e_{" + i + ","  + l + "}");
				}
			}

			// Add constraint c_f:	\sum_{i \in V} \sum_{k \in R_i} v_{ik} y_{ikl} <= S_{l}, j \in P, l \in D_j
			for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
			{
				expr = new GRBLinExpr();
				for(int i = 1; i <= N; i++)			
				{				
					int VM_type_index = getVMTypeIndex(i);
					for(int k = 1; k <= VM_storage_num[VM_type_index]; k++)
						expr.addTerm(VM_storage_size[VM_type_index], y[i][k][l]);
				}
				model.addConstr(expr, GRB.LESS_EQUAL, PM_storage_size[PM_type_index], "c_f_{" + l + "}");
			}

			// Update model

			model.update();

			// Optimize model

			model.optimize();

			// Print objective
			System.out.println("Obj: " + model.get(GRB.DoubleAttr.ObjVal));

			// Print solution
			for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
				System.out.println(z[l].get(GRB.StringAttr.VarName) + " " + z[l].get(GRB.DoubleAttr.X));
			
			int used_disk_num = (int) model.get(GRB.DoubleAttr.ObjVal);
			String returnline = "";
			returnline = returnline + used_disk_num;
			for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
			{
				returnline = returnline + ",";
				long used_disk_size = 0;
				for(int i = 1; i <= N; i++)			
				{
					int VM_type_index = getVMTypeIndex(i);
					for(int k = 1; k <= VM_storage_num[VM_type_index]; k++)
						used_disk_size += y[i][k][l].get(GRB.DoubleAttr.X) * VM_storage_size[VM_type_index];
				}
				returnline = returnline + used_disk_size;
			}
			// Dispose of model and environment

			model.dispose();
			env.dispose();

			return returnline;
		} catch (GRBException e) {
			System.err.println("Error code: " + e.getErrorCode() + ". " + e.getMessage());
			System.exit(0);
		}
		return "";
	}
/*
	public static void main(String[] args)
	{
		if(args.length < 1)
		{
			System.out.println("Usage: Configuration 1~15.");
			System.exit(1);
		}
		
		int PM_type_index = Integer.parseInt(args[0]);
		if(PM_type_index < 1 || PM_type_index > 15)
		{
			System.out.println("Usage: Configuration 1~15.");
			System.exit(1);
		}

		Configuration conf = new Configuration("conf/VM.conf", "conf/PM.conf");
		conf.pre_compute(PM_type_index);
		conf.enumerate_configurations(PM_type_index, "conf/VM_PM_" + args[0] + "_configurations.conf");
	}
*/
	public static void main(String[] args)
	{
		if(args.length < 1)
		{
			System.out.println("Usage: Configuration 1~15.");
			System.exit(1);
		}
		
		int PM_type_index = Integer.parseInt(args[0]);
		if(PM_type_index < 1 || PM_type_index > 15)
		{
			System.out.println("Usage: Configuration 1~15.");
			System.exit(1);
		}

		Configuration conf = new Configuration("conf/VM.conf", "conf/PM.conf");
		conf.compute_resource_assignment("conf/VM_PM_" + args[0] + "_configurations.conf", "conf/VM_PM_" + args[0] + "_resource.conf", Integer.parseInt(args[0]));
	}
}
