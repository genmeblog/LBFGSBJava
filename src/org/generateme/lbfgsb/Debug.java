package org.generateme.lbfgsb;

import java.util.Arrays;
import java.util.Locale;

public class Debug {
	public static boolean DEBUG = false;
	// array and matrix cell formatting
	private static final String cellfmt = "%12s";
	// number formatting
	private static final String numfmt = "%.6f";
	// locale used
	private static final Locale locale = Locale.US;
	// repeat character
	private static final int repeat = 10;

	private static String repeatChar(char c, int n) {
		char[] buff = new char[n];
		Arrays.fill(buff, c);
		return new String(buff);
	}

	public static void debug(char c, String s) {
		String b = repeatChar(c, repeat);
		System.out.println(b + " " + s + " " + b);
	}

	public static void debug(String s) {
		System.out.println(s);
	}

	public static void debug(double[] a) {
		debug(null, a);
	}

	public static void debug(String s, double[] a) {
		if (s != null)
			System.out.print(s);

		System.out.print("[");
		for (double v : a) {
			System.out.format(cellfmt, String.format(locale, numfmt, v));
		}
		System.out.println("]");
	}

	public static void debug(Matrix m) {
		debug(null, m);
	}

	public static void debug(String s, Matrix m) {
		String shift;

		if (s != null) {
			System.out.print(s + "[");
			shift = repeatChar(' ', s.length() + 1);
		} else {
			System.out.print("[");
			shift = " ";
		}

		for (int row = 0; row < m.rows; row++) {
			if (row > 0)
				System.out.print(shift);
			for (int col = 0; col < m.cols; col++) {
				System.out.format(cellfmt, String.format(locale, numfmt, m.get(row, col)));
			}
			if (row == m.rows - 1)
				System.out.print("]");
			System.out.println();
		}
	}

}
