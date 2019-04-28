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

public class CloudAssignment
{
	// data members
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

	int N;				// #VMs
	int M;				// #PMs

	int max_K;			// max num of VM storage disks for any VM Type
	int max_L;			// max num of PM storage disks for any PM Type

	// constructors
	public CloudAssignment(String conf_VM, String conf_PM, String VM_file, String PM_file, boolean bFractional)
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
		for(int i = 1; i <= num_PM_type; i++)
		{
			PM_total += num_PM[i];
			if(PM_total >= PM_index)
			{
				PM_Type_index = i;
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

	public void solve(String output_file)
	{
		//this is an aternative form follows Dr. Xia's original formulation and has additional contrains (8) and (9)
		try {
			GRBEnv	 env	= new GRBEnv("CloudAssignment.log");
			GRBModel  model = new GRBModel(env);

			// Create variables
			// x_{ij}: x[i][j], i = 1 to N, j = 1 to M
			GRBVar [][] x = new GRBVar[N + 1][M + 1];
			
			for(int i = 1; i <= N; i++)			
				for(int j = 1; j <=M; j++)
					x[i][j] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "x_{" + i + "," + j +"}");

			// y_{ikjl}: y[i][k][j][l], i = 1 to N, k 
			GRBVar [][][][] y = new GRBVar[N + 1][max_K + 1][M + 1][max_L + 1];
			for(int i = 1; i <= N; i++)			
			{
				int VM_type_index = getVMTypeIndex(i);
				for(int k = 1; k <= VM_storage_num[VM_type_index]; k++)
				{
					for(int j = 1; j <= M; j++)
					{
						int PM_type_index = getPMTypeIndex(j);
						for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
							y[i][k][j][l] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "y_{" + i + "," + k + "," + j + "," + l + "}");
					}
				}
			}

			GRBVar [] z = new GRBVar[M + 1];
			// z_j: z[j], j = 1 to M
			for(int j = 1; j <=M; j++)
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

			// Add constraint c_a: z_j <= \sum_{i \in V} x_{ij}, j \in P
			for(int j = 1; j <= M; j++)
			{
					expr = new GRBLinExpr();
					expr.addTerm(1.0, z[j]);
					for(int i = 1; i <= N; i++)
						expr.addTerm(-1.0, x[i][j]); 
					model.addConstr(expr, GRB.LESS_EQUAL, 0.0, "c_a_{" + j + "}");
			}

			// Add constraint c_a': B z_j >= \sum_{i \in V} x_{ij}, j \in P
			for(int j = 1; j <= M; j++)
			{
					expr = new GRBLinExpr();
					expr.addTerm(-N, z[j]);
					for(int i = 1; i <= N; i++)
						expr.addTerm(1.0, x[i][j]); 
					model.addConstr(expr, GRB.LESS_EQUAL, 0.0, "c_a'_{" + j + "}");
			}

			// Add constraint c_b: 	y_{ikjl} <= x_{ij}, i \in V, j \in P, k \in R_i, l \in D_j
			for(int i = 1; i <= N; i++)
			{	
				int VM_type_index = getVMTypeIndex(i);
				for(int k = 1; k <= VM_storage_num[VM_type_index]; k++)
				{
					for(int j = 1; j <= M; j++)
					{
						int PM_type_index = getPMTypeIndex(j);
						for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
						{
							expr = new GRBLinExpr();
							expr.addTerm(1.0, y[i][k][j][l]);
							expr.addTerm(-1.0, x[i][j]); 
							model.addConstr(expr, GRB.LESS_EQUAL, 0.0, "c_b_{" + i + "," + k + "," + j + "," + l + "}");
						}
					}
				}
			}

			// Add constraint c_c: 	\sum_{j \in P} \sum_{l \in D_j} y_{ikjl} = 1, i \in V, k \in R_i
			for(int i = 1; i <= N; i++)			
			{
				int VM_type_index = getVMTypeIndex(i);
				for(int k = 1; k <= VM_storage_num[VM_type_index]; k++)
				{
					expr = new GRBLinExpr();
					for(int j = 1; j <= M; j++)
					{
						int PM_type_index = getPMTypeIndex(j);
						for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
						{	
							expr.addTerm(1.0, y[i][k][j][l]);
						}
					}
					model.addConstr(expr, GRB.EQUAL, 1.0, "c_c_{" + i + "," + k + "}");
				}
			}

			// Add constraint c_d: 	\sum_{j \in P} x_{ij} = 1, i \in V
			for(int i = 1; i <= N; i++)
			{
				expr = new GRBLinExpr();
				for(int j = 1; j <= M; j++)
				{		
					expr.addTerm(1.0, x[i][j]); 

				}
				model.addConstr(expr, GRB.EQUAL, 1.0, "c_d_{" + i + "}");
			}	

			// Add constraint c_e:	\sum_{k \in R_i} y_{ikjl} <= 1, i \in V, j \in P, l \in D_j
			for(int i = 1; i <= N; i++)			
			{
				for(int j = 1; j <= M; j++)
				{
					int PM_type_index = getPMTypeIndex(j);
					for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
					{	
						expr = new GRBLinExpr();
						int VM_type_index = getVMTypeIndex(i);
						for(int k = 1; k <= VM_storage_num[VM_type_index]; k++)
							expr.addTerm(1.0, y[i][k][j][l]);
						model.addConstr(expr, GRB.LESS_EQUAL, 1.0, "c_e_{" + i + "," + j + "," + l + "}");
					}
				}
			}

			// Add constraint c_f:	\sum_{i \in V} \sum_{k \in R_i} v_{ik} y_{ikjl} <= S_{jl}, j \in P, l \in D_j
			for(int j = 1; j <= M; j++)
			{
				int PM_type_index = getPMTypeIndex(j);
				for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
				{
					expr = new GRBLinExpr();
					for(int i = 1; i <= N; i++)			
					{				
						int VM_type_index = getVMTypeIndex(i);
						for(int k = 1; k <= VM_storage_num[VM_type_index]; k++)
							expr.addTerm(VM_storage_size[VM_type_index], y[i][k][j][l]);
					}
					model.addConstr(expr, GRB.LESS_EQUAL, PM_storage_size[PM_type_index], "c_f_{" + j + "," + l + "}");
				}
			}

			// Add constraint c_g:	\sum_{i \in V} \alpha_i x_{ij} <= C_j, j \in P
			GRBLinExpr [] CPU_constr = new GRBLinExpr[M + 1];
			for(int j = 1; j <= M; j++)
			{
				CPU_constr[j] = new GRBLinExpr();
				for(int i = 1; i <= N; i++)			
				{	
					int VM_type_index = getVMTypeIndex(i);
					CPU_constr[j].addTerm(VM_CPU[VM_type_index], x[i][j]);
				}
				int PM_type_index = getPMTypeIndex(j);
				model.addConstr(CPU_constr[j], GRB.LESS_EQUAL, PM_CPU[PM_type_index], "c_g_{" + j + "}");
			}

			// Add constraint c_h:	\sum_{i \in V} \beta_i x_{ij} <= M_j, j \in P
			GRBLinExpr [] memory_constr = new GRBLinExpr[M + 1];
			for(int j = 1; j <= M; j++)
			{
				memory_constr[j] = new GRBLinExpr();
				for(int i = 1; i <= N; i++)			
				{				
					int VM_type_index = getVMTypeIndex(i);
					memory_constr[j].addTerm(VM_memory[VM_type_index], x[i][j]);
				}
				int PM_type_index = getPMTypeIndex(j);
				model.addConstr(memory_constr[j], GRB.LESS_EQUAL, PM_memory[PM_type_index], "c_h_{" + j + "}");	
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
			for(int i = 1; i <= N; i++)			
				for(int j = 1; j <= M; j++)
					System.out.println(x[i][j].get(GRB.StringAttr.VarName)
								 + " " + x[i][j].get(GRB.DoubleAttr.X));

			for(int i = 1; i <= N; i++)			
			{
				int VM_type_index = getVMTypeIndex(i);
				for(int k = 1; k <= VM_storage_num[VM_type_index]; k++)
				{
					for(int j = 1; j <= M; j++)
					{
						int PM_type_index = getPMTypeIndex(j);
						for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
							System.out.println(y[i][k][j][l].get(GRB.StringAttr.VarName)
								 + " " + y[i][k][j][l].get(GRB.DoubleAttr.X));
					}
				}
			}
*/
			// z_j: z[j], j = 1 to M
			for(int j = 1; j <= M; j++)
			{
				int PM_type_index = getPMTypeIndex(j);
				System.out.println(z[j].get(GRB.StringAttr.VarName) + " " + z[j].get(GRB.DoubleAttr.X) + ", type = " + PM_type_index);
			}

			// \sum_{i} x_{ij}: #VMs at PM j
			for(int j = 1; j <= M; j++)
			{
				int num_VM_per_PM = 0;
				for(int i = 1; i <= N; i++)
					num_VM_per_PM += x[i][j].get(GRB.DoubleAttr.X);
				System.out.println("#VMs at PM " + j + " = " + num_VM_per_PM);
			}
					
			// used CPUs at PM j
			for(int j = 1; j <= M; j++)
			{
				int PM_type_index = getPMTypeIndex(j);
				System.out.println("#CPUs at PM " + j + " = " + PM_CPU[PM_type_index] + ", #used CPUs = " + CPU_constr[j].getValue());
			}

			// used memory at PM j
			for(int j = 1; j <= M; j++)
			{
				int PM_type_index = getPMTypeIndex(j);
				System.out.println("memory at PM " + j + " = " + PM_memory[PM_type_index] + ", #used memory = " + memory_constr[j].getValue());
			}

			// #used disks at PM j
			for(int j = 1; j <= M; j++)
			{
				int num_used_disks = 0;
				int PM_type_index = getPMTypeIndex(j);
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
				System.out.println("#disks at PM " + j + " = " + PM_storage_num[PM_type_index] + ", #used disks = " + num_used_disks);
			}

			// used disk size at PM j disk l
			for(int j = 1; j <= M; j++)
			{
				int PM_type_index = getPMTypeIndex(j);
				for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
				{
					double used_disk_size = 0;
					for(int i = 1; i <= N; i++)			
					{				
						int VM_type_index = getVMTypeIndex(i);
						for(int k = 1; k <= VM_storage_num[VM_type_index]; k++)
							used_disk_size += VM_storage_size[VM_type_index] * y[i][k][j][l].get(GRB.DoubleAttr.X);
					}
					System.out.println("disk size at PM " + j + " disk " + l + " = " + PM_storage_size[PM_type_index] + ", #used disk size = " + used_disk_size);
				}
			}

			// used disk utility at PM j
			for(int j = 1; j <= M; j++)
			{
				int PM_type_index = getPMTypeIndex(j);
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
				System.out.println("disk utility at PM " + j + " = " + used_disk_size / disk_size);
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
			System.out.println("Usage: CloudAssignment VM_Set.num PM_Set.num output_model.lp beFractional(0/1).");
			System.exit(0);			
		}

		CloudAssignment solver = null;
		if(args[3].compareTo("0") == 0)
			solver = new CloudAssignment("conf/VM.conf", "conf/PM.conf", args[0], args[1], false);
		else
			solver = new CloudAssignment("conf/VM.conf", "conf/PM.conf", args[0], args[1], true);

		solver.solve(args[2]);
	}
}
