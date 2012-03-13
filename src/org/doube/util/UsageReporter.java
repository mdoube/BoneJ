package org.doube.util;

import ij.Prefs;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.*;
import java.nio.charset.Charset;
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
	private static final String utmt = "utmt=event&";
	private static final String utmul = "utmul=" + getLocaleString() + "&";
	private static final String utmje = "utmje=0&";
	private static final String utmfl = "utmfl=11.1%20r102&";
	private static final String utmr = "utmr=-&";
	private static final String utmp = "utmp=%2Fstats&";

	private static String utme = null;
	private static String utmn = null;
	private static String utms = null;
	private static int session = 0;
	private static String utmcc = null;
	private static long firstTime = 0;
	private static long lastTime = 0;
	private static long thisTime = 0;

	private static Random random;

	private static String utmhid;

	private UsageReporter() {
		random = new Random();
	}

	public static UsageReporter reportOn(Object o) {
		utms = "utms=" + session + "&";
		session++;
		utme = "utme=5(Usage*Plugins*" + o.getClass().getName() + ")&";
		utmn = "utmn=" + random.nextInt(Integer.MAX_VALUE) + "&";
		utmhid = "utmhid=" + random.nextInt(Integer.MAX_VALUE) + "&";
		final long time = System.currentTimeMillis() / 1000;
		if (firstTime == 0) {
			firstTime = time;
			thisTime = time;
		}
		lastTime = thisTime;
		thisTime = time;
		utmcc = getCookieString();
		return INSTANCE;
	}

	private static String getCookieString() {
		int cookie = random.nextInt(Integer.MAX_VALUE);
		int randomValue = random.nextInt(Integer.MAX_VALUE);
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
				+ ".103"
				+ "%3B%2B__utmz%3D"
				+ cookie
				+ "."
				+ thisTime
				+ ".79.42.utmcsr%3Dgoogle%7Cutmccn%3D(organic)%7Cutmcmd%3Dorganic%7Cutmctr%3Dbonej%3B";
		return cc;
	}

	public void send() {
		if (!isAllowed())
			return;
		try {
			URL url = new URL(ga + utmwv + utms + utmn + utmhn + utmt + utme
					+ utmcs + utmsr + utmvp + utmsc + utmul + utmje + utmfl
					+ utmdt + utmhid + utmr + utmp + utmac + utmcc);
//			IJ.log(url.toString());
			URLConnection uc = url.openConnection();
			BufferedReader in = new BufferedReader(new InputStreamReader(
					uc.getInputStream()));
			String inputLine;
			while ((inputLine = in.readLine()) != null)
//				IJ.log(inputLine);
				inputLine.length();
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

	private boolean isAllowed() {
		if (!Prefs.get(ReporterOptions.OPTOUTSET, false))
			new ReporterOptions().run("");
		return Prefs.get(ReporterOptions.OPTOUTKEY, true);
	}
}
