package org.doube.util;

/**
 * MultiThreading copyright 2007 Stephan Preibisch
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 2
 * as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */


/**
 * Multithreader utility class for convenient multithreading of ImageJ plugins
 * 
 * @author Stephan Preibisch
 * @author Michael Doube
 * 
 * @see <p>
 *      <a href="http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD"
 *      >http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/
 *      MultiThreading.java;hb=HEAD</a>
 */
public class Multithreader {
	public static void startTask(Runnable run) {
		Thread[] threads = newThreads();

		for (int ithread = 0; ithread < threads.length; ++ithread)
			threads[ithread] = new Thread(run);
		startAndJoin(threads);
	}

	public static void startTask(Runnable run, int numThreads) {
		Thread[] threads = newThreads(numThreads);

		for (int ithread = 0; ithread < threads.length; ++ithread)
			threads[ithread] = new Thread(run);
		startAndJoin(threads);
	}

	public static Thread[] newThreads() {
		int nthread = Runtime.getRuntime().availableProcessors();
		return new Thread[nthread];
	}

	public static Thread[] newThreads(int numThreads) {
		return new Thread[numThreads];
	}

	public static void startAndJoin(Thread[] threads) {
		for (int ithread = 0; ithread < threads.length; ++ithread) {
			threads[ithread].setPriority(Thread.NORM_PRIORITY);
			threads[ithread].start();
		}

		try {
			for (int ithread = 0; ithread < threads.length; ++ithread)
				threads[ithread].join();
		} catch (InterruptedException ie) {
			throw new RuntimeException(ie);
		}
	}
}
