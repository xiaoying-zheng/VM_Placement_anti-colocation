/* Copyright 2013, Gurobi Optimization, Inc. */

/* This example formulates and solves the following cloud assignment model:

	P1: PMs with enumerable configurations; P2: PMs with unenumerable configurations

	min	\sum_{j \in P} z_j a_j 
	s.t.	z_j <= \sum_{i \in V} x_{ij}, j \in P2
                B z_j >= \sum_{i \in V} x_{ij}, j \in P2
		y_{ikjl} <= x_{ij}, i \in V, j \in P2, k \in R_i, l \in D_j
		\sum_{j \in P2} \sum_{l \in D_j} y_{ikjl} = 1, i \in V, k \in R_i
		\sum_{j \in P2} x_{ij} = 1, i \in V
		\sum_{k \in R_i} y_{ikjl} <= 1, i \in V, j \in P2, l \in D_j
		\sum_{i \in V} \sum_{k \in R_i} v_{ik} y_{ikjl} <= S_{jl}, j \in P2, l \in D_j
		\sum_{i \in V} \alpha_i x_{ij} <= C_j, j \in P2
		\sum_{i \in V} \beta_i x_{ij} <= M_j, j \in P2
		x_{ij}, y_{ikjl}, z_j binary

		x'_{jk} <= z_j, j \in P1, k = 1 to num_configurations_total
		\sum_{k} x'_{jt} <= 1, j \in P1
		\sum_{j \in P1} \sum{t} x'_{jt} w_{PM_type_index, t, VM_type_index} >= 
			m_{VM_type_index} - \sum_{i = VM_type_index} \sum_{j \in P2} x_{ij}, VM_type_index \in V
		x'_{jt}, z_j binary
*/

import java.io.*;
import gurobi.*;

public class CloudAssignment_mix
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
	int [] PM_CPU;			// # PM CPUs
	double [] PM_memory;	 	// # PM memory (GB)
	int [] PM_storage_num;		// # PM storage disk
	int [] PM_storage_size;		// # PM storage disk size (GB)
	int [] PM_purchase_cost;	// # fixed cost (normalized)
	int [] PM_operation_cost;	// # fixed cost (normalized)
  	int [] num_PM;			// number of PMs of each PM types
	boolean [] PM_conf_usage;	// whether use pre-computed configuration; false: not use; true: use

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

	int num_PM_P1 = 0;		//#PMs in set P1
	int num_PM_P2 = 0;		//#PMs in set P2

	// constructors
	public CloudAssignment_mix(String conf_VM, String conf_PM, String VM_file, String PM_file, boolean bFractional)
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
		configurations = new int[num_PM_type + 1][max_num_configurations + 1][num_VM_type + 1];	// configurations: configurations[PM_type_index][conf_index][VM_type_index]: #VMs 

		//Read PM_conf_usage configuration
		try{
			bufferedReader = new BufferedReader(new FileReader("conf/PM_conf_usage.conf"));
			PM_conf_usage = new boolean[num_PM_type + 1];
			for(int j = 1; j <= num_PM_type; j++)
			{
				//usage
				line = bufferedReader.readLine();
				if(Integer.parseInt(line) == 1){
					PM_conf_usage[j] = true;
					num_PM_P1 += num_PM[j];
				}else{
					PM_conf_usage[j] = false;
					num_PM_P2 += num_PM[j];
				}
			}
			bufferedReader.close();
			System.out.println("num_PM_P1 = " + num_PM_P1 + ", num_PM_P2 = " + num_PM_P2);
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
		int VM_type_index = 0;

		int VM_total = 0;		
		for(int i = 1; i <= num_VM_type; i++)
		{
			VM_total += num_VM[i];
			if(VM_total >= VM_index)
			{
				VM_type_index = i;
				break;
			}
		}

		if(0 == VM_type_index)
		{
			System.out.println("Error: wrong VM_index!");
			System.exit(1);
		}

		return VM_type_index;
	}

	// input: PM index; output: PM type index;
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
			System.out.println("Error: getPMTypeIndex(): wrong PM_index = " + PM_index + " !");
			System.exit(1);
		}
		return PM_Type_index;
	}

	// input: PM index; output: PM type index; for those in the set P1
	private int getPMTypeIndex_P1(int PM_index)
	{
		int PM_Type_index = 0;

		int PM_total = 0;		
		for(int j = 1; j <= num_PM_type; j++)
		{
			if(PM_conf_usage[j] == false)
				continue;

			PM_total += num_PM[j];
			if(PM_total >= PM_index)
			{
				PM_Type_index = j;
				break;
			}
		}
		if(0 == PM_Type_index)
		{
			System.out.println("Error: getPMTypeIndex_P1() wrong PM_index = " + PM_index + " !");
			System.exit(1);
		}
		return PM_Type_index;
	}

	// input: PM index; output: PM type index; for those in the set P2
	private int getPMTypeIndex_P2(int PM_index)
	{
		int PM_Type_index = 0;

		int PM_total = 0;		
		for(int j = 1; j <= num_PM_type; j++)
		{
			if(PM_conf_usage[j] == true)
				continue;

			PM_total += num_PM[j];
			if(PM_total >= PM_index)
			{
				PM_Type_index = j;
				break;
			}
		}
		if(0 == PM_Type_index)
		{
			System.out.println("Error: getPMTypeIndex_P2(): wrong PM_index = " + PM_index + " !");
			System.exit(1);
		}
		return PM_Type_index;
	}

	// input: index in the set P1; output: PM index
	private int getPMIndexfromP1Index(int P1_index)
	{
		int PM_Type_index = getPMTypeIndex_P1(P1_index);

		int PM_index = P1_index;

		int j = 1;
		while(j < PM_Type_index)
		{
			if(PM_conf_usage[j] == false)
				PM_index += num_PM[j];				
			j++;
		}
		return PM_index;
	}

	// input: index in the set P2; output: PM index
	private int getPMIndexfromP2Index(int P2_index)
	{
		int PM_Type_index = getPMTypeIndex_P2(P2_index);

		int PM_index = P2_index;

		int j = 1;
		while(j < PM_Type_index)
		{
			if(PM_conf_usage[j] == true)
				PM_index += num_PM[j];				
			j++;
		}
		return PM_index;	
	}

	// input: VM type index; output: lower index of VMs
	private int getVMLowerIndex(int VM_type_index)
	{
		int VM_lower_index = 1;

		int j = 1;
		while(j < VM_type_index)
		{
			VM_lower_index += num_VM[j];				
			j++;
		}
		return VM_lower_index;	
	}

	// input: VM type index; output: upper index of VMs
	private int getVMUpperIndex(int VM_type_index)
	{
		int VM_upper_index = 0;

		int j = 1;
		while(j <= VM_type_index)
		{
			VM_upper_index += num_VM[j];				
			j++;
		}
		return VM_upper_index;	
	}

	private void PM_print(int PM_P1_index, int configuration_index, int option)
	{
		/* option = 1: print #VMs
		   option = 2: print #used CPUs
		   option = 3: print used memory
		   option = 4: print #used disks
		   option = 5: print used disks size per disk

		*/

		int PM_type_index = getPMTypeIndex_P1(PM_P1_index);
		int j = getPMIndexfromP1Index(PM_P1_index);
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
		try {
			GRBEnv	 env	= new GRBEnv("CloudAssignment_mix.log");
			GRBModel  model = new GRBModel(env);

			GRBVar [][] x = null;
			GRBVar [][] x_prime = null;
			GRBVar [][][][] y = null;

			if(num_PM_P2 > 0){
				// Create variables
				// x_{ij}: x[i][j], i = 1 to N, j \in P2
				x = new GRBVar[N + 1][num_PM_P2 + 1];
			
				for(int i = 1; i <= N; i++)			
					for(int j = 1; j <= num_PM_P2; j++)
						x[i][j] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "x_{" + i + "," + j +"}");
			}

			if(num_PM_P1 > 0){
				// x'_{jt}: x'[j][t], j = \in P1, t = 1 to #configurations
				x_prime = new GRBVar[num_PM_P1 + 1][max_num_configurations + 1];
			
				for(int j = 1; j <= num_PM_P1; j++)			
				{
					int PM_type_index = getPMTypeIndex_P1(j);
					for(int t = 1; t <= num_configurations[PM_type_index]; t++)
						x_prime[j][t] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "x_prime{" + j + "," + t +"}");
				}

			}

			if(num_PM_P2 > 0){
				// y_{ikjl}: y[i][k][j][l], i = 1 to N, k, j \in P_2
				y = new GRBVar[N + 1][max_K + 1][num_PM_P2 + 1][max_L + 1];
				for(int i = 1; i <= N; i++)			
				{
					int VM_type_index = getVMTypeIndex(i);
					for(int k = 1; k <= VM_storage_num[VM_type_index]; k++)
					{
						for(int j = 1; j <= num_PM_P2; j++)
						{
							int PM_type_index = getPMTypeIndex_P2(j);
							for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
								y[i][k][j][l] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "y_{" + i + "," + k + "," + j + "," + l + "}");
						}
					}
				}
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

			GRBLinExpr [] CPU_constr = null;
			GRBLinExpr [] memory_constr = null;
			if(num_PM_P2 > 0){
				// Add constraint c_a: z_j <= \sum_{i \in V} x_{ij}, j \in P_2
				for(int j = 1; j <= num_PM_P2; j++)
				{
					expr = new GRBLinExpr();
					int index_P = getPMIndexfromP2Index(j);
					expr.addTerm(1.0, z[index_P]);
					for(int i = 1; i <= N; i++)
						expr.addTerm(-1.0, x[i][j]); 
					model.addConstr(expr, GRB.LESS_EQUAL, 0.0, "c_a_{" + j + "}");
				}
	
				// Add constraint c_b: B z_j >= \sum_{i \in V} x_{ij}, j \in P_2
				for(int j = 1; j <= num_PM_P2; j++)
				{
					expr = new GRBLinExpr();
					int index_P = getPMIndexfromP2Index(j);
					expr.addTerm(-N, z[index_P]);
					for(int i = 1; i <= N; i++)
						expr.addTerm(1.0, x[i][j]); 
					model.addConstr(expr, GRB.LESS_EQUAL, 0.0, "c_b_{" + j + "}");
				}

				// Add constraint c_c: 	y_{ikjl} <= x_{ij}, i \in V, j \in P_2, k \in R_i, l \in D_j
				for(int i = 1; i <= N; i++)
				{	
					int VM_type_index = getVMTypeIndex(i);
					for(int k = 1; k <= VM_storage_num[VM_type_index]; k++)
					{
						for(int j = 1; j <= num_PM_P2; j++)
						{
							int PM_type_index = getPMTypeIndex_P2(j);
							for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
							{
								expr = new GRBLinExpr();
								expr.addTerm(1.0, y[i][k][j][l]);
								expr.addTerm(-1.0, x[i][j]); 
								model.addConstr(expr, GRB.LESS_EQUAL, 0.0, "c_c_{" + i + "," + k + "," + j + "," + l + "}");
							}
						}
					}
				}

				// Add constraint c_d: 	\sum_{j \in P_2} \sum_{l \in D_j} y_{ikjl} <= \sum_{j \in P_2} x_{ij}, i \in V, k \in R_i
				for(int i = 1; i <= N; i++)			
				{
					int VM_type_index = getVMTypeIndex(i);
					for(int k = 1; k <= VM_storage_num[VM_type_index]; k++)
					{
						expr = new GRBLinExpr();
						for(int j = 1; j <= num_PM_P2; j++)
						{
							int PM_type_index = getPMTypeIndex_P2(j);
							for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
								expr.addTerm(1.0, y[i][k][j][l]);

							expr.addTerm(-1.0, x[i][j]);
						}

						model.addConstr(expr, GRB.EQUAL, 0.0, "c_d_{" + i + "," + k + "}");
					}
				}

				// Add constraint c_e: 	\sum_{j \in P_2} x_{ij} <= 1, i \in V
				for(int i = 1; i <= N; i++)
				{
					expr = new GRBLinExpr();
					for(int j = 1; j <= num_PM_P2; j++)
						expr.addTerm(1.0, x[i][j]); 
					model.addConstr(expr, GRB.LESS_EQUAL, 1.0, "c_e_{" + i + "}");
				}	

				// Add constraint c_f:	\sum_{k \in R_i} y_{ikjl} <= 1, i \in V, j \in P_2, l \in D_j
				for(int i = 1; i <= N; i++)			
				{
					for(int j = 1; j <= num_PM_P2; j++)
					{
						int PM_type_index = getPMTypeIndex_P2(j);
						for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
						{	
							expr = new GRBLinExpr();
							int VM_type_index = getVMTypeIndex(i);
							for(int k = 1; k <= VM_storage_num[VM_type_index]; k++)
								expr.addTerm(1.0, y[i][k][j][l]);
							model.addConstr(expr, GRB.LESS_EQUAL, 1.0, "c_f_{" + i + "," + j + "," + l + "}");
						}
					}
				}

				// Add constraint c_g:	\sum_{i \in V} \sum_{k \in R_i} v_{ik} y_{ikjl} <= S_{jl}, j \in P_2, l \in D_j
				for(int j = 1; j <= num_PM_P2; j++)
				{
					int PM_type_index = getPMTypeIndex_P2(j);
					for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
					{
						expr = new GRBLinExpr();
						for(int i = 1; i <= N; i++)			
						{				
							int VM_type_index = getVMTypeIndex(i);
							for(int k = 1; k <= VM_storage_num[VM_type_index]; k++)
								expr.addTerm(VM_storage_size[VM_type_index], y[i][k][j][l]);
						}
						model.addConstr(expr, GRB.LESS_EQUAL, PM_storage_size[PM_type_index], "c_g_{" + j + "," + l + "}");
					}
				}

				// Add constraint c_h:	\sum_{i \in V} \alpha_i x_{ij} <= C_j, j \in P_2
				CPU_constr = new GRBLinExpr[num_PM_P2 + 1];
				for(int j = 1; j <= num_PM_P2; j++)
				{
					CPU_constr[j] = new GRBLinExpr();
					for(int i = 1; i <= N; i++)			
					{	
						int VM_type_index = getVMTypeIndex(i);
						CPU_constr[j].addTerm(VM_CPU[VM_type_index], x[i][j]);
					}
					int PM_type_index = getPMTypeIndex_P2(j);
					model.addConstr(CPU_constr[j], GRB.LESS_EQUAL, PM_CPU[PM_type_index], "c_h_{" + j + "}");
				}

				// Add constraint c_i:	\sum_{i \in V} \beta_i x_{ij} <= M_j, j \in P_2
				memory_constr = new GRBLinExpr[num_PM_P2 + 1];
				for(int j = 1; j <= num_PM_P2; j++)
				{
					memory_constr[j] = new GRBLinExpr();
					for(int i = 1; i <= N; i++)			
					{				
						int VM_type_index = getVMTypeIndex(i);
						memory_constr[j].addTerm(VM_memory[VM_type_index], x[i][j]);
					}
					int PM_type_index = getPMTypeIndex_P2(j);
					model.addConstr(memory_constr[j], GRB.LESS_EQUAL, PM_memory[PM_type_index], "c_i_{" + j + "}");	
				}
			}
/************************************************************************************************/

			if(num_PM_P1 > 0){
				// Add constraint c_j: z_j <= \sum_{j \in P_1} x'_{jt}, t = 1 to ...
				for(int j = 1; j <= num_PM_P1; j++)
				{
					int PM_type_index = getPMTypeIndex_P1(j);
					expr = new GRBLinExpr();
					int index_P = getPMIndexfromP1Index(j);
					expr.addTerm(1.0, z[index_P]);
					for(int t = 1; t <= num_configurations[PM_type_index]; t++)
						expr.addTerm(-1.0, x_prime[j][t]); 
					model.addConstr(expr, GRB.LESS_EQUAL, 0.0, "c_j_{" + j + "}");
				}

				// Add constraint c_k: 	B z_j >= \sum_{j \in P_1} x'_{jt}, t = 1 to ...; B = max_num_configurations
				for(int j = 1; j <= num_PM_P1; j++)
				{
					int PM_type_index = getPMTypeIndex_P1(j);
					expr = new GRBLinExpr();
					int index_P = getPMIndexfromP1Index(j);
					expr.addTerm(-max_num_configurations, z[index_P]);
					for(int t = 1; t <= num_configurations[PM_type_index]; t++)
						expr.addTerm(1.0, x_prime[j][t]); 
					model.addConstr(expr, GRB.LESS_EQUAL, 0.0, "c_k_{" + j + "}");
				}

				// Add constraint c_l: 	\sum_{t} x'_{jt} <= 1, \forall j = 1 to num_PM_P1
				for(int j = 1; j <= num_PM_P1; j++)
				{	
					int PM_type_index = getPMTypeIndex_P1(j);
					expr = new GRBLinExpr();
					for(int t = 1; t <= num_configurations[PM_type_index]; t++)
						expr.addTerm(1.0, x_prime[j][t]);
					model.addConstr(expr, GRB.LESS_EQUAL, 1.0, "c_l_{" + j + "}");
				}

				// Add constraint c_m: 	\sum_{j \in P_1} \sum{t} x'_{jt} w_{PM_type_index, t, VM_type_index} >= m_{VM_type_index}  - \sum_{i = VM_type_index} \sum_{j \in P_2} x_{ij}, \forall VM_type_index \in V
				for(int VM_type_index = 1; VM_type_index <= num_VM_type; VM_type_index++)			
				{
					expr = new GRBLinExpr();
					for(int j = 1; j <= num_PM_P1; j++)
					{
						int PM_type_index = getPMTypeIndex_P1(j);
						for(int t = 1; t <= num_configurations[PM_type_index]; t++)
							expr.addTerm(configurations[PM_type_index][t][VM_type_index], x_prime[j][t]);
					}
	
					int VM_lower_index = getVMLowerIndex(VM_type_index);
					int VM_upper_index = getVMUpperIndex(VM_type_index);
					for(int i = VM_lower_index; i <= VM_upper_index; i++)
					{
						for(int j = 1; j <= num_PM_P2; j++)
							expr.addTerm(1.0, x[i][j]);				
					}

					model.addConstr(expr, GRB.GREATER_EQUAL, num_VM[VM_type_index], "c_m_{" + VM_type_index + "}");
				}
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
			// z_j: z[j], j = 1 to M
			for(int j = 1; j <= M; j++)
				System.out.println(z[j].get(GRB.StringAttr.VarName)
								 + " " + z[j].get(GRB.DoubleAttr.X));

			// print #VMs
			for(int j = 1; j <= num_PM_P1; j++)			
			{
				int PM_type_index = getPMTypeIndex_P1(j);
				boolean flag = false;
				for(int t = 1; t <= num_configurations[PM_type_index]; t++)
				{
					if(x_prime[j][t].get(GRB.DoubleAttr.X) == 1.0)
					{
						flag = true;
						PM_print(j, t, 1);
						break;
					}
				}
				if(flag == false)
					System.out.println("#VMs at PM " + getPMIndexfromP1Index(j) + " = 0");
			}

			// \sum_{i} x_{ij}: #VMs at PM j
			for(int j = 1; j <= num_PM_P2; j++)
			{
				int num_VM_per_PM = 0;
				for(int i = 1; i <= N; i++)
					num_VM_per_PM += x[i][j].get(GRB.DoubleAttr.X);
				System.out.println("#VMs at PM " + getPMIndexfromP2Index(j) + " = " + num_VM_per_PM);
			}

			// used CPUs at PM j
			for(int j = 1; j <= num_PM_P1; j++)
			{
				int PM_type_index = getPMTypeIndex_P1(j);
				boolean flag = false;
				for(int t = 1; t <= num_configurations[PM_type_index]; t++)
				{
					if(x_prime[j][t].get(GRB.DoubleAttr.X) == 1)
					{
						flag = true;
						PM_print(j, t, 2);
						break;
					}
				}
				if(flag == false)
					System.out.println("#CPUs at PM " + getPMIndexfromP1Index(j) + " = " + PM_CPU[PM_type_index] + ", #used CPUs = 0");
			}

			// used CPUs at PM j
			for(int j = 1; j <= num_PM_P2; j++)
			{
				int PM_type_index = getPMTypeIndex_P2(j);
				System.out.println("#CPUs at PM " + getPMIndexfromP2Index(j) + " = " + PM_CPU[PM_type_index] + ", #used CPUs = " + CPU_constr[j].getValue());
			}

			// used memory at PM j
			for(int j = 1; j <= num_PM_P1; j++)
			{
				int PM_type_index = getPMTypeIndex_P1(j);
				boolean flag = false;
				for(int t = 1; t <= num_configurations[PM_type_index]; t++)
				{
					if(x_prime[j][t].get(GRB.DoubleAttr.X) == 1)
					{
						flag = true;
						PM_print(j, t, 3);
						break;
					}
				}
				if(flag == false)
					System.out.println("memory at PM " + getPMIndexfromP1Index(j) + " = " + PM_memory[PM_type_index] + ", #used memory = 0");
			}

			// used memory at PM j
			for(int j = 1; j <= num_PM_P2; j++)
			{
				int PM_type_index = getPMTypeIndex_P2(j);
				System.out.println("memory at PM " + getPMIndexfromP2Index(j) + " = " + PM_memory[PM_type_index] + ", #used memory = " + memory_constr[j].getValue());
			}

			// #used disks at PM j
			for(int j = 1; j <= num_PM_P1; j++)
			{
				int PM_type_index = getPMTypeIndex_P1(j);
				boolean flag = false;
				for(int t = 1; t <= num_configurations[PM_type_index]; t++)
				{
					if(x_prime[j][t].get(GRB.DoubleAttr.X) == 1)
					{
						flag = true;
						PM_print(j, t, 4);
						break;
					}
				}
				if(flag == false)
					System.out.println("#disks at PM " + getPMIndexfromP1Index(j) + " = " + PM_storage_num[PM_type_index] + ", #used disks = 0");
			}

			// #used disks at PM j
			for(int j = 1; j <= num_PM_P2; j++)
			{
				int num_used_disks = 0;
				int PM_type_index = getPMTypeIndex_P2(j);
				for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
				{
					boolean flag = false;
					for(int i = 1; i <= N; i++)			
					{
						int VM_type_index = getVMTypeIndex(i);
						for(int k = 1; k <= VM_storage_num[VM_type_index]; k++)
						{
							if(y[i][k][j][l].get(GRB.DoubleAttr.X) == 1.0) 
							{
								num_used_disks += 1;
								flag = true;
								break;
							}
						}
						if(flag == true)
							break;
					}
				}
				System.out.println("#disks at PM " + getPMIndexfromP2Index(j) + " = " + PM_storage_num[PM_type_index] + ", #used disks = " + num_used_disks);
			}


			// used disk size at PM j disk l
			for(int j = 1; j <= num_PM_P1; j++)
			{
				int PM_type_index = getPMTypeIndex_P1(j);
				boolean flag = false;
				for(int t = 1; t <= num_configurations[PM_type_index]; t++)
				{
					if(x_prime[j][t].get(GRB.DoubleAttr.X) == 1)
					{
						flag = true;
						PM_print(j, t, 5);
						break;
					}
				}
				if(flag == false)
				{
					for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
						System.out.println("disk size at PM " + getPMIndexfromP1Index(j) + " disk " + l + " = " + PM_storage_size[PM_type_index] + ", #used disk size = 0");
				}

			}

			// used disk size at PM j disk l
			for(int j = 1; j <= num_PM_P2; j++)
			{
				int PM_type_index = getPMTypeIndex_P2(j);
				for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
				{
					double used_disk_size = 0;
					for(int i = 1; i <= N; i++)			
					{				
						int VM_type_index = getVMTypeIndex(i);
						for(int k = 1; k <= VM_storage_num[VM_type_index]; k++)
							used_disk_size += VM_storage_size[VM_type_index] * y[i][k][j][l].get(GRB.DoubleAttr.X);
					}
					System.out.println("disk size at PM " + getPMIndexfromP2Index(j) + " disk " + l + " = " + PM_storage_size[PM_type_index] + ", #used disk size = " + used_disk_size);
				}
			}

			// used disk utility at PM j disk l
			for(int j = 1; j <= num_PM_P1; j++)
			{
				int PM_type_index = getPMTypeIndex_P1(j);
				boolean flag = false;
				for(int t = 1; t <= num_configurations[PM_type_index]; t++)
				{
					if(x_prime[j][t].get(GRB.DoubleAttr.X) == 1)
					{
						flag = true;
						PM_print(j, t, 6);
						break;
					}
				}
				if(flag == false)
					System.out.println("disk utility at PM " + getPMIndexfromP1Index(j) + " = 0");
			}

			// used disk utility at PM j
			for(int j = 1; j <= num_PM_P2; j++)
			{
				int PM_type_index = getPMTypeIndex_P2(j);
				double disk_size = 0;
				double used_disk_size = 0;
				for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
				{
					for(int i = 1; i <= N; i++)			
					{				
						int VM_type_index = getVMTypeIndex(i);
						for(int k = 1; k <= VM_storage_num[VM_type_index]; k++)
							used_disk_size += VM_storage_size[VM_type_index] * y[i][k][j][l].get(GRB.DoubleAttr.X);
					}
					disk_size += PM_storage_size[PM_type_index];
				}
				System.out.println("disk utility at PM " + getPMIndexfromP2Index(j) + " = " + used_disk_size / disk_size);
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
			System.out.println("Usage: CloudAssignment_mix VM_Set.num PM_Set.num output_model.lp beFractional(0/1).");
			System.exit(0);			
		}

		CloudAssignment_mix solver = null;
		if(args[3].compareTo("0") == 0)
			solver = new CloudAssignment_mix("conf/VM.conf", "conf/PM.conf", args[0], args[1], false);
		else
			solver = new CloudAssignment_mix("conf/VM.conf", "conf/PM.conf", args[0], args[1], true);

		solver.solve(args[2]);
	}
}
