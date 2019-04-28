all: CloudAssignment.class CloudAssignment_form1.class CloudAssignment_mix.class Configuration.class RandomAssignment.class RandomAssignment_restricted.class
CloudAssignment.class: CloudAssignment.java
	javac -classpath /opt/gurobi652/linux64/lib/gurobi.jar:. CloudAssignment.java
CloudAssignment_form1.class: CloudAssignment_form1.java
	javac -classpath /opt/gurobi652/linux64/lib/gurobi.jar:. CloudAssignment_form1.java
CloudAssignment_mix.class: CloudAssignment_mix.java
	javac -classpath /opt/gurobi652/linux64/lib/gurobi.jar:. CloudAssignment_mix.java
Configuration.class: Configuration.java
	javac -classpath /opt/gurobi652/linux64/lib/gurobi.jar:. Configuration.java
RandomAssignment.class: RandomAssignment.java
	javac RandomAssignment.java
RandomAssignment_restricted.class: RandomAssignment_restricted.java
	javac RandomAssignment_restricted.java
clean:
	rm *.class *.log *~
