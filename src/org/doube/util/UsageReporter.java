package org.doube.util;

import ij.IJ;//only for debug logging

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.*;
import java.nio.charset.Charset;
import java.util.Date;
import java.util.Locale;
import java.util.Random;

public class UsageReporter {
	public static final UsageReporter INSTANCE = new UsageReporter();

	private static final String ga = "http://www.google-analytics.com/__utm.gif?";
	private static final String utmwv = "utmwv=5.2.5&";
	private static final String utmhn = "utmhn=bonej.org&";
	private static final String utmcs = "utmcs=" + Charset.defaultCharset()
			+ "&";
	private static final String utmsr = "utmsr=1280x800&";
	private static final String utmvp = "utmvp=1280x690&";
	private static final String utmsc = "utmsc=24-bit&";
	private static final String utmac = "utmac=UA-366405-8&";
	private static final String utmdt = "utmdt=bonej.org%20Usage%20Statistics&";
	private static final String utmt = "utmt=page&"; // better, event
	private static final String utmul = "utmul=" + getLocaleString() + "&";
	private static final String utmje = "utmje=0&";
	private static final String utmfl = "utmfl=11.1%20r102&";
	// private static final String utmhid = "utmhid=1811992293&"; random adsense
	private static final String utmr = "utmr=0&";

	private static String utmc = null;
	private static String utmn = null;
	private static String utmp = null;
	private static int utms = -1;
	private static String utmcc = null;
	private static long firstTime = 0;
	private static long lastTime = 0;
	private static long thisTime = 0;

	private static Random random;
	private static Date date;

	private UsageReporter() {
		random = new Random();
		date = new Date();
	}

	public static UsageReporter getInstance(Object o) {
		utms++;
		// handle repeated calls from the same user session
		if (utmc == null) {
			utmc = "utmcn=1&";
		} else {
			utmc = "utmcr=1&";
		}
		// set a new utmn per request to avoid caching of the gif
		utmn = "utmn="
				+ Integer
						.toString((int) (Math.floor(Math.random() * 999999999)))
				+ "&";

		utmp = "utmp=%2Fstats?caller=" + o.getClass().getName() + "&";
		utmcc = getCookieString();
		send();
		return INSTANCE;
	}

	private static String getCookieString() {
		int cookie = random.nextInt(Integer.MAX_VALUE);
		int randomValue = random.nextInt(2147483647);
		lastTime = thisTime;
		thisTime = date.getTime();
		if (firstTime == 0) {
			firstTime = thisTime;
			lastTime = thisTime;
		}
		String cc = "utmcc=__utma%3D"
				+ cookie
				+ "."
				+ randomValue
				+ "."

				+ firstTime
				+ "."
				+ lastTime
				+ "."
				+ thisTime
				+ ".3"
				+ "%3B%2B__utmb%3D"
				+ cookie
				+

				"%3B%2B__utmc%3D"
				+ cookie
				+

				"%3B%2b__utmz%3D"
				+ cookie
				+ "."
				+ thisTime
				+ ".3.2.utmccn%3D(organic)%7Cutmcsr%3Dgoogle%7Cutmctr%3Dbonej%2Busage%2Blogging%7Cutmcmd%3Dorganic%3B%2B";
		return cc;
	}

	private static void send() {
		try {
			URL url = new URL(ga + utmwv + utms + utmn + utmhn + utmcs + utmsr
					+ utmvp + utmsc + utmul + utmje + utmfl + utmr + utmp
					+ utmac + utmdt + utmc + utmt + utmcc);
			IJ.log(url.toString());
			URLConnection uc = url.openConnection();
			BufferedReader in = new BufferedReader(new InputStreamReader(
					uc.getInputStream()));
			String inputLine;
			while ((inputLine = in.readLine()) != null)
				IJ.log(inputLine);
			in.close();
		} catch (MalformedURLException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static String getLocaleString() {
		String locale = Locale.getDefault().toString();
		locale = locale.replace("_", "-");
		locale = locale.toLowerCase();
		return locale;
	}
}
