/* Copyright 2013, Gurobi Optimization, Inc. */

/* This example formulates and solves the following cloud assignment model:

	min	\sum_{j \in P} z_j a_j 
	s.t.	x_{jk} <= z_j, j \in P, k = 1 ...
		\sum_{k} x_{jk} <= 1, j \in P
		\sum_{j \in P} \sum{k} x_{jk} w_{PM_type_index, k, VM_type_index} >= m_{VM_type_index}, VM_type_index \in V
		x_{jk}, z_j binary
*/

import java.io.*;
import gurobi.*;

public class CloudAssignment_form1
{
	// data members
	final int max_num_configurations = 4247;	// max number of configurations for any PM type. pre-computed.
	int num_VM_type;		// number of VM types
	String [] VM_type;  		// VM types
	int [] VM_CPU;			// # VM CPUs
	double [] VM_memory; 		// # VM memory (GB)
	int [] VM_storage_num;		// # VM storage disk
	int [] VM_storage_size;		// # VM storage disk size (GB)
	int [] num_VM;			// number of VMs of each VM types

	int num_PM_type;		// number of PM types
	String [] PM_type;  		// PM types
	int [] PM_CPU;			// # VM CPUs
	double [] PM_memory;	 	// # VM memory (GB)
	int [] PM_storage_num;		// # VM storage disk
	int [] PM_storage_size;		// # VM storage disk size (GB)
	int [] PM_purchase_cost;	// # fixed cost (normalized)
	int [] PM_operation_cost;	// # fixed cost (normalized)
  	int [] num_PM;			// number of PMs of each PM types

	int [] num_configurations;	// number of configurations of each PM type
	int [][][] configurations;	// configurations: configurations[PM_type_index][conf_index][VM_type_index]: #VMs 
	int [][] allocated_VMs;		// allocated_VMs: allocated_VMs[PM_type_index][conf_index]
	int [][] used_CPUs;		// used_CPUs: used_CPUs[PM_type_index][conf_index]
	double [][] used_memory;	// used_memory: used_memory[PM_type_index][conf_index]    
	int [][] used_disk_num;		// used_disk_num: configurations[PM_type_index][conf_index]
	int [][][] used_disk_size;	// used_disk_size: configurations[PM_type_index][conf_index][disk_index] 

	int N;				// #VMs
	int M;				// #PMs

	int max_K;			// max num of VM storage disks for any VM Type
	int max_L;			// max num of PM storage disks for any PM Type

	// constructors
	public CloudAssignment_form1(String conf_VM, String conf_PM, String VM_file, String PM_file, boolean bFractional)
	{
		String line = null;
		N = 0;
		M = 0;
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
			num_VM = new int[num_VM_type + 1];

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

			//#VMs, read from VM.num file. if bFractional == true, num = N * fraction
			bufferedReader = new BufferedReader(new FileReader(VM_file));
			
			//#VMs
			line = bufferedReader.readLine();

			//#total: N
			line = bufferedReader.readLine();
			int total_VM = Integer.parseInt(line);
			
			for(int i = 1; i <= num_VM_type; i++)
			{
				line = bufferedReader.readLine();
				if(bFractional)
					num_VM[i] = (int)(total_VM * Double.parseDouble(line));
				else
					num_VM[i] = Integer.parseInt(line);
				N += num_VM[i];
			}

			bufferedReader.close();

			if(N != total_VM)
				System.out.println("Warning: N = " + N + ", but the input value = " + total_VM);
		}catch(IOException e){
			System.out.println("Error: " + e.toString());
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
			num_PM = new int[num_PM_type + 1];
		
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

			//#PMs, read from PM.num file. if bFractional == true, num = M * fraction
			bufferedReader = new BufferedReader(new FileReader(PM_file));
			
			//#PMs
			line = bufferedReader.readLine();

			//#total: M
			line = bufferedReader.readLine();
			int total_PM = Integer.parseInt(line);
			
			for(int i = 1; i <= num_PM_type; i++)
			{
				line = bufferedReader.readLine();
				if(bFractional)
					num_PM[i] = (int)(total_PM * Double.parseDouble(line));
				else
					num_PM[i] = Integer.parseInt(line);
				M += num_PM[i];
			}

			bufferedReader.close();

			if(M != total_PM)
				System.out.println("Warning: M = " + M + ", but the input value = " + total_PM);
		}catch(IOException e){
			System.out.println("Error: " + e.toString());
		}

		num_configurations = new int[num_PM_type + 1];		// number of configurations of each PM type
		// configurations: configurations[PM_type_index][conf_index][VM_type_index]: #VMs 
		configurations = new int[num_PM_type + 1][max_num_configurations + 1][num_VM_type + 1];
		// allocated_VMs: allocated_VMs[PM_type_index][conf_index]
		allocated_VMs = new int[num_PM_type + 1][max_num_configurations + 1];
		// used_CPUs: used_CPUs[PM_type_index][conf_index]
		used_CPUs = new int[num_PM_type + 1][max_num_configurations + 1];
		// used_memory: used_memory[PM_type_index][conf_index]
		used_memory = new double[num_PM_type + 1][max_num_configurations + 1];
		// used_disk_num: used_disk_num[PM_type_index][conf_index]: #used_disk
		used_disk_num = new int[num_PM_type + 1][max_num_configurations + 1];
		// used_disk_size: used_disk_size[PM_type_index][conf_index][disk_index]: used_disk_size
		used_disk_size = new int[num_PM_type + 1][max_num_configurations + 1][PM_storage_num[num_PM_type] + 1];	

		//Read PM-VM configuration
		try{
			String configuration_file;

			for(int j = 1; j <= num_PM_type; j++)
			{
				configuration_file = "conf/VM_PM_" + j + "_configurations.conf";
				bufferedReader = new BufferedReader(new FileReader(configuration_file));
				System.out.println("conf/VM_PM_" + j + "_configurations.conf");

				//#configurations
				line = bufferedReader.readLine();
				num_configurations[j] = Integer.parseInt(line);

				for(int index = 1; index <= num_configurations[j]; index++)
				{
					line = bufferedReader.readLine();
					int pre_pos = 0;
					int pos = 0;
					for(int i = 1; i <= num_VM_type - 1; i++)
					{
						pos = line.indexOf(",", pre_pos);
						configurations[j][index][i] = Integer.parseInt(line.substring(pre_pos, pos));
						pre_pos = pos + 1;
					}
					configurations[j][index][num_VM_type] = Integer.parseInt(line.substring(pre_pos));
				}

				bufferedReader.close();
				
				// Read PM resource assignment
				configuration_file = "conf/VM_PM_" + j + "_resource.conf";
				bufferedReader = new BufferedReader(new FileReader(configuration_file));
				System.out.println("conf/VM_PM_" + j + "_resource.conf");

				for(int index = 1; index <= num_configurations[j]; index++)
				{
					// #VMs,#CPUs,memory,disk num,disk size
					line = bufferedReader.readLine();
					int pre_pos = 0;
					int pos = line.indexOf(",", pre_pos);
					allocated_VMs[j][index] = Integer.parseInt(line.substring(pre_pos, pos));

					pre_pos = pos + 1;
					pos = line.indexOf(",", pre_pos);
					used_CPUs[j][index] = Integer.parseInt(line.substring(pre_pos, pos));

					pre_pos = pos + 1;
					pos = line.indexOf(",", pre_pos);
					used_memory[j][index] = Double.parseDouble(line.substring(pre_pos, pos));

					pre_pos = pos + 1;
					pos = line.indexOf(",", pre_pos);
					used_disk_num[j][index] = Integer.parseInt(line.substring(pre_pos, pos));

					pre_pos = pos + 1;
					for(int l = 1; l <= PM_storage_num[j] - 1; l++)
					{
						
						pos = line.indexOf(",", pre_pos);
						used_disk_size[j][index][l] = Integer.parseInt(line.substring(pre_pos, pos));
						pre_pos = pos + 1;
					}
					used_disk_size[j][index][PM_storage_num[j]] = Integer.parseInt(line.substring(pre_pos));
				}

				bufferedReader.close();
			}
		}catch(IOException e){
			System.out.println("Error: " + e.toString());
		}
	}

	// methods

	// input: VM index; output: VM type index
	private int getVMTypeIndex(int VM_index)
	{
		int VM_Type_index = 0;

		int VM_total = 0;		
		for(int i = 1; i <= num_VM_type; i++)
		{
			VM_total += num_VM[i];
			if(VM_total >= VM_index)
			{
				VM_Type_index = i;
				break;
			}
		}

		if(0 == VM_Type_index)
		{
			System.out.println("Error: wrong VM_index!");
			System.exit(1);
		}

		return VM_Type_index;
	}

	// input: PM index; output: PM type index
	private int getPMTypeIndex(int PM_index)
	{
		int PM_Type_index = 0;

		int PM_total = 0;		
		for(int j = 1; j <= num_PM_type; j++)
		{
			PM_total += num_PM[j];
			if(PM_total >= PM_index)
			{
				PM_Type_index = j;
				break;
			}
		}
		if(0 == PM_Type_index)
		{
			System.out.println("Error: wrong PM_index!");
			System.exit(1);
		}
		return PM_Type_index;
	}
 
	private void PM_print(int PM_index, int configuration_index, int option)
	{
		/* option = 1: print #VMs
		   option = 2: print #used CPUs
		   option = 3: print used memory
		   option = 4: print #used disks
		   option = 5: print used disks size per disk

		*/

		int j = PM_index;
		int PM_type_index = getPMTypeIndex(PM_index);
		int k = configuration_index;
		if(1 == option){
			// #VMs at PM j
			System.out.println("#VMs at PM " + j + " = " + allocated_VMs[PM_type_index][k]);
		}else if(2 == option){
			// used CPUs at PM j
			System.out.println("#CPUs at PM " + j + " = " + PM_CPU[PM_type_index] + ", #used CPUs = " + used_CPUs[PM_type_index][k]);
		}else if(3 == option){
			// used memory at PM j
			System.out.println("memory at PM " + j + " = " + PM_memory[PM_type_index] + ", #used memory = " + used_memory[PM_type_index][k]);
		}else if(4 == option){
			// #used disks at PM j
			System.out.println("#disks at PM " + j + " = " + PM_storage_num[PM_type_index] + ", #used disks = " + used_disk_num[PM_type_index][k]);
		}else if(5 == option){
			// used disk size at PM j
			for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
				System.out.println("disk size at PM " + j + " disk " + l + " = " + PM_storage_size[PM_type_index] + ", #used disk size = " + used_disk_size[PM_type_index][k][l]);
		}else if(6 == option){
			// used disk utility at PM j
			double disk_size = 0;
			double disk_size_used = 0;
			for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
			{
				disk_size += PM_storage_size[PM_type_index];
				disk_size_used += used_disk_size[PM_type_index][k][l];
			}
			System.out.println("disk utility at PM " + j + " = " + disk_size_used / disk_size);
		}

	}

	public void solve(String output_file)
	{
/* This example formulates and solves the following cloud assignment model:

	min	\sum_{j \in P} z_j a_j 
	s.t.	z_j <= \sum_{j \in P} x_{jk}, k = 1 to ...
		B z_j >= \sum_{j \in P} x_{jk}, k = 1 to ...; B = max_num_configurations
		\sum_{k} x_{jk} <= 1, j \in P
		\sum_{j \in P} \sum{k} x_{jk} w_{PM_type_index, k, VM_type_index} >= m_{VM_type_index}, VM_type_index \in V
		x_{jk}, z_j binary
*/


		try {
			GRBEnv	 env	= new GRBEnv("CloudAssignment_form1.log");
			GRBModel  model = new GRBModel(env);

			// Create variables
			// x_{jk}: x[j][k], j = 1 to M, k = 1 to max_num_configurations
			GRBVar [][] x = new GRBVar[M + 1][max_num_configurations + 1];
			
			for(int j = 1; j <= M; j++)			
			{
				int PM_type_index = getPMTypeIndex(j);
				for(int k = 1; k <= num_configurations[PM_type_index]; k++)
					x[j][k] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "x_{" + j + "," + k +"}");
			}

			GRBVar [] z = new GRBVar[M + 1];
			// z_j: z[j], j = 1 to M
			for(int j = 1; j <= M; j++)
				z[j] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "z_" + j);

			// Integrate new variables
			model.update();

			// Set objective: min	\sum_{j \in P} z_j a_j
			GRBLinExpr expr = new GRBLinExpr();
			for(int j = 1; j <= M; j++)
			{
				int PM_type_index = getPMTypeIndex(j);

				expr.addTerm(PM_operation_cost[PM_type_index], z[j]);
			}
			model.setObjective(expr, GRB.MINIMIZE);

			// Add constraint c_a: z_j <= \sum_{j \in P} x_{jk}, k = 1 to ...
			for(int j = 1; j <= M; j++)
			{
				int PM_type_index = getPMTypeIndex(j);
				expr = new GRBLinExpr();
				expr.addTerm(1.0, z[j]);
				for(int k = 1; k <= num_configurations[PM_type_index]; k++)
					expr.addTerm(-1.0, x[j][k]); 
				model.addConstr(expr, GRB.LESS_EQUAL, 0.0, "c_a_{" + j + "}");
			}

			// Add constraint c_b: 	B z_j >= \sum_{j \in P} x_{jk}, k = 1 to ...; B = max_num_configurations
			for(int j = 1; j <= M; j++)
			{
				int PM_type_index = getPMTypeIndex(j);
				expr = new GRBLinExpr();
				expr.addTerm(-max_num_configurations, z[j]);
				for(int k = 1; k <= num_configurations[PM_type_index]; k++)
					expr.addTerm(1.0, x[j][k]); 
				model.addConstr(expr, GRB.LESS_EQUAL, 0.0, "c_b_{" + j + "}");
			}

			// Add constraint c_c: 	\sum_{k} x_{jk} <= 1, \forall j = 1 to M
			for(int j = 1; j <= M; j++)
			{	
				int PM_type_index = getPMTypeIndex(j);
				expr = new GRBLinExpr();
				for(int k = 1; k <= num_configurations[PM_type_index]; k++)
					expr.addTerm(1.0, x[j][k]);
				model.addConstr(expr, GRB.LESS_EQUAL, 1.0, "c_b_{" + j + "}");
			}

			// Add constraint c_d: 	\sum_{j \in P} \sum{k} x_{jk} w_{PM_type_index, k, VM_type_index} >= m_{VM_type_index}, \forall VM_type_index \in V
			for(int VM_type_index = 1; VM_type_index <= num_VM_type; VM_type_index++)			
			{
				expr = new GRBLinExpr();
				for(int j = 1; j <= M; j++)
				{
					int PM_type_index = getPMTypeIndex(j);
					for(int k = 1; k <= num_configurations[PM_type_index]; k++)
						expr.addTerm(configurations[PM_type_index][k][VM_type_index], x[j][k]);
				}
				model.addConstr(expr, GRB.GREATER_EQUAL, num_VM[VM_type_index], "c_c_{" + VM_type_index + "}");
			}

			// Update model

			model.update();

			// write the model to a file
			model.write(output_file);

			// Optimize model

			model.optimize();

			// Print objective
			System.out.println("Obj: " + model.get(GRB.DoubleAttr.ObjVal));

			// Print solution
/*
			for(int j = 1; j <= M; j++)			
			{
				int PM_type_index = getPMTypeIndex(j);
				for(int k = 1; k <= num_configurations[PM_type_index]; k++)
					System.out.println(x[j][k].get(GRB.StringAttr.VarName)
								 + " " + x[j][k].get(GRB.DoubleAttr.X));
			}
*/


			// z_j: z[j], j = 1 to M
			for(int j = 1; j <= M; j++)
			{
				int PM_type_index = getPMTypeIndex(j);
				System.out.println(z[j].get(GRB.StringAttr.VarName) + " " + z[j].get(GRB.DoubleAttr.X) + ", type = " + PM_type_index);
			}
			
			// #VMs at PM j
			for(int j = 1; j <= M; j++)			
			{
				int PM_type_index = getPMTypeIndex(j);
				boolean flag = false;
				for(int k = 1; k <= num_configurations[PM_type_index]; k++)
				{
					if(x[j][k].get(GRB.DoubleAttr.X) == 1)
					{
						flag = true;
						PM_print(j, k, 1);
						break;
					}
				}
				if(flag == false)
					System.out.println("#VMs at PM " + j + " = 0");
			}

			// used CPUs at PM j
			for(int j = 1; j <= M; j++)
			{
				int PM_type_index = getPMTypeIndex(j);
				boolean flag = false;
				for(int k = 1; k <= num_configurations[PM_type_index]; k++)
				{
					if(x[j][k].get(GRB.DoubleAttr.X) == 1)
					{
						flag = true;
						PM_print(j, k, 2);
						break;
					}
				}
				if(flag == false)
					System.out.println("#CPUs at PM " + j + " = " + PM_CPU[PM_type_index] + ", #used CPUs = 0");
			}

			// used memory at PM j
			for(int j = 1; j <= M; j++)
			{
				int PM_type_index = getPMTypeIndex(j);
				boolean flag = false;
				for(int k = 1; k <= num_configurations[PM_type_index]; k++)
				{
					if(x[j][k].get(GRB.DoubleAttr.X) == 1)
					{
						flag = true;
						PM_print(j, k, 3);
						break;
					}
				}
				if(flag == false)
					System.out.println("memory at PM " + j + " = " + PM_memory[PM_type_index] + ", #used memory = 0");
			}

			// #used disks at PM j
			for(int j = 1; j <= M; j++)
			{
				int PM_type_index = getPMTypeIndex(j);
				boolean flag = false;
				for(int k = 1; k <= num_configurations[PM_type_index]; k++)
				{
					if(x[j][k].get(GRB.DoubleAttr.X) == 1)
					{
						flag = true;
						PM_print(j, k, 4);
						break;
					}
				}
				if(flag == false)
					System.out.println("#disks at PM " + j + " = " + PM_storage_num[PM_type_index] + ", #used disks = 0");
			}

			// used disk size at PM j disk l
			for(int j = 1; j <= M; j++)
			{
				int PM_type_index = getPMTypeIndex(j);
				boolean flag = false;
				for(int k = 1; k <= num_configurations[PM_type_index]; k++)
				{
					if(x[j][k].get(GRB.DoubleAttr.X) == 1)
					{
						flag = true;
						PM_print(j, k, 5);
						break;
					}
				}
				if(flag == false)
				{
					for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
						System.out.println("disk size at PM " + j + " disk " + l + " = " + PM_storage_size[PM_type_index] + ", #used disk size = 0");
				}
			}

			// used disk utility at PM j disk l
			for(int j = 1; j <= M; j++)
			{
				int PM_type_index = getPMTypeIndex(j);
				boolean flag = false;
				for(int k = 1; k <= num_configurations[PM_type_index]; k++)
				{
					if(x[j][k].get(GRB.DoubleAttr.X) == 1)
					{
						flag = true;
						PM_print(j, k, 6);
						break;
					}
				}
				if(flag == false)
					System.out.println("disk utility at PM " + j + " = 0");
			}

			// Dispose of model and environment

			model.dispose();
			env.dispose();

		 } catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
								 e.getMessage());
		 }
	}


	public static void main(String[] args)
	{
		if(args.length != 4)
		{
			System.out.println("Usage: CloudAssignment_form1 VM_Set.num PM_Set.num output_model.lp beFractional(0/1).");
			System.exit(0);			
		}

		CloudAssignment_form1 solver = null;
		if(args[3].compareTo("0") == 0)
			solver = new CloudAssignment_form1("conf/VM.conf", "conf/PM.conf", args[0], args[1], false);
		else
			solver = new CloudAssignment_form1("conf/VM.conf", "conf/PM.conf", args[0], args[1], true);

		solver.solve(args[2]);
	}
}
