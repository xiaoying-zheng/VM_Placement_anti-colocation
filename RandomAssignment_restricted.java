/* Copyright 2013, Gurobi Optimization, Inc. */

/* This example formulates and solves the following cloud assignment model with a random solution

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



permute the VM_list
for each VM in VM_list do
	scan used_PM_list to find a PM that can accommodate the VM
	if such a PM is found then
		assign the VM to the PM
	else
		scan unused_PM_list to find a PM that can accommodate the VM
		if such a PM is found then
			assign the VM to the PM; move the PM to	used_PM_list
		else
			exit		# the problem is infeasible
		end if
	end if
end for
*/

import java.io.*;
import java.util.*;

public class RandomAssignment_restricted
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

	int [] assignment;		// assignment[i] = j; VM i is assigned to PM j

	boolean [] PM_ON;		// PM_ON[j] == true if PM j is turned on
	int [] PM_CPU_balance;		// CPU balance of PM[j]
	double [] PM_memory_balance;	// memory balance of PM[j]
	int [][] PM_storage_balance;	// storage balance of PM[j]
	
	int max_K;			// max num of VM storage disks for any VM Type
	int max_L;			// max num of PM storage disks for any PM Type

	// constructors
	public RandomAssignment_restricted(String conf_VM, String conf_PM, String VM_file, String PM_file, boolean bFractional)
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

			assignment = new int[N + 1];
			
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

			PM_ON = new boolean[M + 1];
			PM_CPU_balance = new int[M + 1];
			PM_memory_balance = new double[M + 1];
			PM_storage_balance = new int[M + 1][max_L + 1];

			for(int i = 1; i <= M; i++)
			{
				int PM_type_index = getPMTypeIndex(i);

				PM_ON[i] = false;
				PM_CPU_balance[i] = PM_CPU[PM_type_index];
				PM_memory_balance[i] = PM_memory[PM_type_index];

				int storage_num = PM_storage_num[PM_type_index];				

				for(int j = 1; j <= max_L; j++)
				{
					if(j <= storage_num)
						PM_storage_balance[i][j] = PM_storage_size[PM_type_index];
					else
						PM_storage_balance[i][j] = 0;
				}
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

	// check if PM j can accomodate VM type i
	public boolean checkAccomodation(int PM_index, int VM_type_index)
	{

		int PM_type_index = getPMTypeIndex(PM_index);

		if(PM_type_index == 10)
		{
			if(VM_type_index <= 2)
				return false;
		}else if(PM_type_index == 11 || PM_type_index == 12){
		
			if(VM_type_index <= 3)
				return false;
			if(VM_type_index == 5 || VM_type_index == 6 || VM_type_index == 10 || VM_type_index == 11 || VM_type_index == 15)
				return false;
		}else if(PM_type_index >= 13){
			if(VM_type_index <= 7 || VM_type_index == 10 || VM_type_index == 11 || VM_type_index == 12 || VM_type_index == 15 || VM_type_index == 16)
				return false;
		}

		if(PM_CPU_balance[PM_index] < VM_CPU[VM_type_index])
			return false;
		if(PM_memory_balance[PM_index] < VM_memory[VM_type_index])
			return false;
		

		int qualified_disk_num = 0;
		for(int l = 1; l <= max_L; l++)
			if(PM_storage_balance[PM_index][l] >= VM_storage_size[VM_type_index])
				qualified_disk_num++;

		if(qualified_disk_num < VM_storage_num[VM_type_index])
			return false;

		return true;
	}

	public void assignVM(int PM_index, int VM_index, int VM_type_index)
	{
		assignment[VM_index] = PM_index;
		PM_ON[PM_index] = true;
		
		PM_CPU_balance[PM_index] -= VM_CPU[VM_type_index];

		PM_memory_balance[PM_index] -= VM_memory[VM_type_index];
		
		int l = 1;
		for(int k = 1; k <= VM_storage_num[VM_type_index]; k++)
		{
			while(l <= max_L)
			{
				if(PM_storage_balance[PM_index][l] >= VM_storage_size[VM_type_index])
				{
					PM_storage_balance[PM_index][l] -= VM_storage_size[VM_type_index];
					l++;
					break;
				}
				l++;
			}
		}	
	}

	public double solve(String output_file)
	{

		for(int i = 1; i <= M; i++)
		{
			int PM_type_index = getPMTypeIndex(i);

			PM_ON[i] = false;
			PM_CPU_balance[i] = PM_CPU[PM_type_index];
			PM_memory_balance[i] = PM_memory[PM_type_index];

			int storage_num = PM_storage_num[PM_type_index];				

			for(int j = 1; j <= max_L; j++)
			{
				if(j <= storage_num)
					PM_storage_balance[i][j] = PM_storage_size[PM_type_index];
				else
					PM_storage_balance[i][j] = 0;
			}
		}

		List<Integer> list = new ArrayList<Integer>();
		for(int i = 1; i <= N; i++)
			list.add(i);
		java.util.Collections.shuffle(list);

		double obj = 0;

		for(int i = 1; i <= N; i++)
		{
			int VM_index = list.get(i - 1);

			int VM_type_index = getVMTypeIndex(VM_index);

			int selected_PM_index = -1;
			
			// scan the used PMs	
			for(int j = 1; j <= M; j++)
			{
				if(true == PM_ON[j])
					if(checkAccomodation(j, VM_type_index))
					{					
						selected_PM_index = j;
						break;
					}
			}

			if(-1 != selected_PM_index)
			{
				assignVM(selected_PM_index, VM_index, VM_type_index);
				continue;
			}

			// scan the unused PMs
			for(int j = 1; j <= M; j++)
			{
				if(false == PM_ON[j])
					if(checkAccomodation(j, VM_type_index))
					{				
						selected_PM_index = j;
						break;
					}
			}

			if(-1 != selected_PM_index)
			{
				obj += PM_operation_cost[getPMTypeIndex(selected_PM_index)];
				assignVM(selected_PM_index, VM_index, VM_type_index);
				continue;
			}
			System.err.println("Unable to assign VM = " + VM_index + ", VM_type_index = " + VM_type_index);
		}			

		// Print objective
		System.out.println("Obj: " + obj);
/*
		// Print solution
		System.out.println("VM-PM assignment:");

		for(int i = 1; i <= N; i++)			
			System.out.println("VM = " + i + ", PM = " + assignment[i]);

		System.out.println("PM turn on/off:");

		int PM_ON_num = 0;
		for(int j = 1; j <= M; j++)
		{
			int PM_type_index = getPMTypeIndex(j);
			System.out.println("PM = " + j + ", PM_ON[j] = " + PM_ON[j] + ", type = " + PM_type_index);
			if(PM_ON[j] == true)			
				PM_ON_num++;
		}

		System.out.println("#PM ON = " + PM_ON_num);
			
		for(int j = 1; j <= M; j++)
		{
			if(PM_ON[j] == false)
				continue;

			int PM_type_index = getPMTypeIndex(j);

			System.out.print("PM index = " + j + ", type = " + PM_type_index);

			// used CPUs percent at PM j
			System.out.print(", CPU = " + (double) (PM_CPU[PM_type_index] - PM_CPU_balance[j]) / PM_CPU[PM_type_index]);

			// used memory at PM j
			System.out.print(", memory = " + (PM_memory[PM_type_index] - PM_memory_balance[j]) / PM_memory[PM_type_index]);

			// #used disks at PM j
			int num_used_disks = 0;
			for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
			{
				if(PM_storage_balance[j][l] < PM_storage_size[PM_type_index])
					num_used_disks++;

			}
			System.out.print(", #disks = " + (double) num_used_disks / PM_storage_num[PM_type_index]);

			// max used disk utility at PM j
			double max_disk_size_percent = 0;
			for(int l = 1; l <= PM_storage_num[PM_type_index]; l++)
			{
				double used_disk_size_percent = (double)(PM_storage_size[PM_type_index] - PM_storage_balance[j][l]) / PM_storage_size[PM_type_index];
				if(max_disk_size_percent < used_disk_size_percent)
					max_disk_size_percent = used_disk_size_percent;
			}
			System.out.println(", disk utility = " + max_disk_size_percent);
		}
*/
		return obj;
	}

	public static void main(String[] args)
	{
		if(args.length != 4)
		{
			System.out.println("Usage: RandomAssignment_restricted VM_Set.num PM_Set.num output_model.lp beFractional(0/1).");
			System.exit(0);			
		}

		RandomAssignment_restricted solver = null;
		if(args[3].compareTo("0") == 0)
			solver = new RandomAssignment_restricted("conf/VM.conf", "conf/PM.conf", args[0], args[1], false);
		else
			solver = new RandomAssignment_restricted("conf/VM.conf", "conf/PM.conf", args[0], args[1], true);

//		solver.solve(args[2]);

		double obj = 0;
		for(int loop = 1; loop <= 1000; loop++)
			obj += solver.solve(args[2]);

		System.out.println("Average obj = " + obj / 1000);
	}
}
