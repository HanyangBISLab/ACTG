package Environments;

public class Error {
	
	public static final int PARAM_ERROR = 0;
	public static final int NO_SUCH_A_FILE = 1;
	public static final int IS_NOT_A_DIRECTORY = 2;
	public static final int IS_EMPTY = 3;
	public static final int FORMAT_ERROR = 4;
	public static final int FACTAL_ERROR = 5;
	
	public static void exitError(int status, String msg){
		System.out.println("*-----------------------ERROR MESSAGE-------------------*");
		
		if(status == PARAM_ERROR){
			System.out.println("The params file is invalid: " + msg);
		}else if(status == NO_SUCH_A_FILE){
			System.out.println("There is no such a file: " + msg);
		}else if(status == IS_NOT_A_DIRECTORY){
			System.out.println("It is not a directory: " + msg);
		}else if(status == IS_EMPTY){
			System.out.println("The directory is empty: " + msg);
		}else if(status == FORMAT_ERROR){
			System.out.println("The format is invalid: " + msg);
		}else if(status == FACTAL_ERROR){
			System.out.println("The factal error is occurred: " + msg);
		}
		
		System.out.println("Please see the README file.");
		System.exit(-1);
	}
}