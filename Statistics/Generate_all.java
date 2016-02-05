
public class Generate_all {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		int size = 15;
		Generate_n gen = new Generate_n();
		String[] s = new String[1];
		for(int i = 0; i < 55+1; i++){
			s[0] = Integer.toString(15 + i);
			gen.main(s);
			System.out.println((15+i) + " done.");
		}
	}
}
