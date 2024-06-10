package main.java.util;

import main.java.FDR.FDREval;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class Cloger {

	private static Cloger instance = null;
	private long start_time;

	public Logger logger;


	private Cloger(){
		start_time = System.currentTimeMillis();
		logger = LogManager.getLogger(FDREval.class.getName());
	}

	public static Cloger getInstance() {
		if (instance == null) {
			instance = new Cloger();
		}
		return instance;
	}
}


