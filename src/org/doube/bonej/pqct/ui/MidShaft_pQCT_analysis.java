/*
	This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

	N.B.  the above text was copied from http://www.gnu.org/licenses/gpl.html
	unmodified. I have not attached a copy of the GNU license to the source...

    Copyright (C) 2011 Timo Rantalainen
*/

/*
Written by Timo Rantalainen tjrantal@gmail.com March 2011

Contributions from other people taken from the internet
Cubic Spline equations:					http://www.physics.arizona.edu/~restrepo/475A/Notes/sourcea/node35.html
linear equation group solver:			http://math.nist.gov/javanumerics/jama/ 			Used for solving the spline equations
erode algorithm (modified from dilate):	http://ostermiller.org/dilate_and_erode.html

The file structure of Stratec files was volunteered from Stratec in response to my query

*/
package ui;
import javax.swing.*;		//GUI using swing
import java.awt.event.*; 	//Events and Actionlisteners
import java.io.*;				//File IO
import java.lang.Math;
import java.awt.*;
import java.awt.geom.Line2D;
import javax.swing.event.*;
import javax.swing.border.*;
import java.util.Vector;
import java.util.Enumeration;
import java.io.*;
import javax.sound.sampled.*;
import java.awt.font.*;
import java.text.*;
import java.awt.image.*;
import java.awt.image.DataBuffer;
import Jama.*; //For solving linear equation groups, i.e. for  calculating the spline http://math.nist.gov/javanumerics/jama/
import Analysis.*;	//Analysis stuff..
import SelectRoi.*;	//ROI selection..
import ImageReading.*;	//image data (includes DICOMReader, which seems to work for uncompressed explicit and implicit VR-files.  lossles JPEG not implemented)
import DrawImage.*;		//Drawing and saving images
/*implements AL enables using events. Operating system will implement...*/
/*extends = inherit, object can inherit only one object*/
public class MidShaft_pQCT_analysis extends JPanel implements ActionListener, ItemListener {	
	JLabel bmd;
	JLabel area;
	JLabel marrowCenter;
	JButton setFolderIn;
	JButton setSaveFolder;
	JButton setSaveImageFolder;
	JButton analyze;
	JButton erode;
	JButton save;
	
	JCheckBox manualRoiDetermination;
	JCheckBox distalBone;
	JCheckBox erodeOn;
	JCheckBox marrowAna;
	JCheckBox cortAna;
	JCheckBox distributionAna;
	JCheckBox softAna;
	JCheckBox ukkSpecial;
	JCheckBox dicomAna;
	JCheckBox imSave;
	JTextField softThresholdInput;
	JTextField boneThresholdInput;
	JTextField marrowThresholdInput;
	JTextField areaThresholdInput;
	JTextField BMDThresholdInput;
	JTextField filterInput;
	JTextField softFilterInput;
	JTextField tNimi;
	DrawImage showFigure;
	JTextField tScaleCoeff;
	JTextField tScaleConst;

	int sectorWidth;
	String currentPath;
	String savePath;
	String imageSavePath;
	File fileIn;
	MouseListener coordinateWindowPressed;
	Vector<Double> roiX;
	Vector<Double> roiY;
	
	double softThreshold;// = Double.parseDouble(softThresholdInput.getText());
	double boneThreshold;// = Double.parseDouble(boneThresholdInput.getText());
	double marrowThreshold;// = Double.parseDouble(marrowThresholdInput.getText());
	double areaThreshold; //= Double.parseDouble(areaThresholdInput.getText());
	double BMDThreshold; //= Double.parseDouble(BMDThresholdInput.getText());

	
	int filterSize;
	int softFilterSize;
	boolean mRoiDet;	//Manual ROI determination
	boolean dBoneSite;	//Is distal bone site going to be analyzed
	boolean eOn;	//Erode for distal bone sites (should not be used ever..)
	boolean cOn;	//Cortical analysis
	boolean mOn;	//Marrow analysis
	boolean dOn;	//Density distribution analysis
	boolean stOn;	//Soft tissue analysis
	boolean ukkOn;	//UKK special to remove the plastic sleeve ring from XCT3000 images measured at UKK institute
	boolean dicomOn;	//Whether dicom files are to be analyzed
	boolean imOn;		//Whether we want images saved or not
	DecimalFormat oneDec;
	
	ImageAndAnalysisDetails details;
	public static void doAndShowGUI(){
	/*
	        JFrame frame = new JFrame("MouseMotionEventDemo");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        
        //Create and set up the content pane.
        JComponent newContentPane = new MouseMotionEventDemo();
        newContentPane.setOpaque(true); //content panes must be opaque
        frame.setContentPane(newContentPane);
        
        //Display the window.
        frame.pack();
        frame.setVisible(true);
	*/
	
		MidShaft_pQCT_analysis midshaft_pQCT_analysis = new MidShaft_pQCT_analysis();
		JFrame f = new JFrame("pQCT Analysis");
		f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		JComponent newContentPane = new MidShaft_pQCT_analysis();
		newContentPane.setOpaque(true); //content panes must be opaque
        f.setContentPane(newContentPane);

        f.pack();
        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
        int w = 750;	//Width of frame on screen
        int h = 700;	//Height of frame on screen
        f.setLocation(screenSize.width/2 - w/2, screenSize.height/2 - h/2);
        f.setSize(w, h);
        f.setVisible(true);
		
	}
	
	public MidShaft_pQCT_analysis(){
		JPanel buttons = new JPanel();
		setFolderIn= new JButton("pQCT path");
		setFolderIn.setMnemonic(KeyEvent.VK_A);
		setFolderIn.setActionCommand("setFolderIn");
		setFolderIn.addActionListener(this);
		setFolderIn.setToolTipText("Press to set pQCT image path.");
		buttons.add(setFolderIn);
		setSaveFolder= new JButton("save path");
		setSaveFolder.setMnemonic(KeyEvent.VK_B);
		setSaveFolder.setActionCommand("setSaveFolder");
		setSaveFolder.addActionListener(this);
		setSaveFolder.setToolTipText("Press to set result file path.");
		buttons.add(setSaveFolder);
		setSaveImageFolder= new JButton("image save path");
		setSaveImageFolder.setMnemonic(KeyEvent.VK_C);
		setSaveImageFolder.setActionCommand("setSaveImageFolder");
		setSaveImageFolder.addActionListener(this);
		setSaveImageFolder.setToolTipText("Press to set result file path.");
		buttons.add(setSaveImageFolder);
		manualRoiDetermination = new JCheckBox("ManualROI",false);
        manualRoiDetermination.setMnemonic(KeyEvent.VK_D);
		manualRoiDetermination.addItemListener(this);
		buttons.add(manualRoiDetermination);
		mRoiDet = false;
		distalBone = new JCheckBox("DistalBone",false);
        distalBone.setMnemonic(KeyEvent.VK_E);
		distalBone.addItemListener(this);
		buttons.add(distalBone);
		dBoneSite = false;
		

		analyze = new JButton("Analyze");
		analyze.setMnemonic(KeyEvent.VK_F);
		analyze.setActionCommand("analyze");
		analyze.addActionListener(this);
		analyze.setToolTipText("Press to analyze.");
		buttons.add(analyze);
		add(buttons);
		/*
		erode = new JButton("Erode");
		erode.setMnemonic(KeyEvent.VK_C);
		erode.setActionCommand("erode");
		erode.addActionListener(this);
		erode.setToolTipText("Show erode.");
		buttons.add(erode);
		*/
		/*
		save = new JButton("Save");
		save.setMnemonic(KeyEvent.VK_D);
		save.setActionCommand("save");
		save.addActionListener(this);
		save.setToolTipText("Show erode.");
		buttons.add(save);
		*/
		JPanel selections = new JPanel();
		
		erodeOn		= new JCheckBox("erodeOn",false);
        erodeOn.setMnemonic(KeyEvent.VK_G);
		erodeOn.addItemListener(this);
		selections.add(erodeOn);
		eOn =false;
		
		marrowAna = new JCheckBox("MarrowAna",false);
        marrowAna.setMnemonic(KeyEvent.VK_H);
		marrowAna.addItemListener(this);
		selections.add(marrowAna);
		mOn = false;
		distributionAna= new JCheckBox("DistributionAna",true);
        distributionAna.setMnemonic(KeyEvent.VK_I);
		distributionAna.addItemListener(this);
		selections.add(distributionAna);
		dOn = true;
		
		softAna = new JCheckBox("SoftTissueAna",false);
        softAna.setMnemonic(KeyEvent.VK_J);
		softAna.addItemListener(this);
		selections.add(softAna);
		stOn = false;
		
		cortAna= new JCheckBox("CortAna",true);
        cortAna.setMnemonic(KeyEvent.VK_K);
		cortAna.addItemListener(this);
		selections.add(cortAna);
		cOn =true;
		
		ukkSpecial= new JCheckBox("ukkSpecial",false);
        ukkSpecial.setMnemonic(KeyEvent.VK_L);
		ukkSpecial.addItemListener(this);
		selections.add(ukkSpecial);
		ukkOn =false;
		
		dicomAna= new JCheckBox("DICOM",false);
        dicomAna.setMnemonic(KeyEvent.VK_M);
		dicomAna.addItemListener(this);
		selections.add(dicomAna);
		dicomOn = false;
		
		imSave= new JCheckBox("Visual",true);
        imSave.setMnemonic(KeyEvent.VK_N);
		imSave.addItemListener(this);
		selections.add(imSave);
		imOn = true;
		
		selections.setOpaque(true); //content panes must be opaque
		add(selections);
		JPanel textSelections = new JPanel();
		
		
		/*
		JTextField softThresholdInput;
		JTextField boneThresholdInput;
		JTextField marrowThresholdInput;
		JTextField areaThresholdInput;
		JTextField BMDThresholdInput;
		JTextField filterInput;
		JTextField softFilterInput;
		*/
		
		softThresholdInput = new JTextField("300");
		//softThresholdInput = new JTextField("1400");
		textSelections.add(softThresholdInput);
		softThreshold = Double.parseDouble(softThresholdInput.getText());
		boneThresholdInput = new JTextField("550");
		//boneThresholdInput = new JTextField("1600");
		textSelections.add(boneThresholdInput);
		boneThreshold = Double.parseDouble(boneThresholdInput.getText());
		marrowThresholdInput = new JTextField("300");
		//marrowThresholdInput = new JTextField("1400");
		textSelections.add(marrowThresholdInput);
		marrowThreshold = Double.parseDouble(marrowThresholdInput.getText());
		areaThresholdInput = new JTextField("550");
		//areaThresholdInput = new JTextField("1600");
		textSelections.add(areaThresholdInput);
		areaThreshold = Double.parseDouble(areaThresholdInput.getText());
		BMDThresholdInput = new JTextField("690");
		//BMDThresholdInput = new JTextField("1800");
		textSelections.add(BMDThresholdInput);
		BMDThreshold = Double.parseDouble(BMDThresholdInput.getText());
		filterInput = new JTextField("3");
		textSelections.add(filterInput);
		filterSize = Integer.parseInt(filterInput.getText());
		softFilterInput = new JTextField("7");
		textSelections.add(softFilterInput);
		softFilterSize = Integer.parseInt(softFilterInput.getText());
		tNimi = new JTextField("marrowDist.xls");
		textSelections.add(tNimi);
//	tScaleCoeff = 1.495;		//pQCT 3000
//	tScaleConst = -341.0;			//pQCT 3000
//	tScaleCoeff = 1.724;		//pQCT 2000
//	tScaleConst = -322.0;			//pQCT 2000
		//tScaleCoeff = new JTextField("1.495");
		tScaleCoeff = new JTextField("1.724");
		//tScaleCoeff = new JTextField("1.0");
		textSelections.add(tScaleCoeff);
		
		//tScaleConst = new JTextField("-341.0");
		tScaleConst = new JTextField("-322.0");
		//tScaleConst = new JTextField("0.0");
		textSelections.add(tScaleConst);
		textSelections.setOpaque(true); //content panes must be opaque
		add(textSelections);
		//ADD textfields..

		showFigure = new DrawImage();
		showFigure.setBackground(new Color(0, 0, 0));
		showFigure.setPreferredSize(new Dimension(500,500));
		showFigure.setOpaque(true);
		
		add(showFigure);
		showFigure.addMouseListener(showFigure);
		
		currentPath = new String();
		savePath = new String();
		imageSavePath = new String();
		
		/*Preset path*/
		//File vali = new File("/home/timo/Deakin/LOOK/LOOK_pQCT_Scans");
		/*
		File vali = new File("C:/Oma/Deakin/LOOK/LOOK_pQCT_Scans");
		currentPath  =vali.getAbsolutePath();
		savePath = vali.getAbsolutePath();
		imageSavePath = vali.getAbsolutePath();
		*/
		/*CURRENT PATH*/
		
		currentPath  =System.getProperty("user.dir");
		savePath = System.getProperty("user.dir");
		imageSavePath = System.getProperty("user.dir");
		
		oneDec = new DecimalFormat("0");
	}
	
	public void actionPerformed(ActionEvent e) {
		if ("setFolderIn".equals(e.getActionCommand())){
			JFileChooser chooser = new JFileChooser(currentPath);
			chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
			int returnVal = chooser.showOpenDialog(this);
			String pathIn = new String();
			fileIn = new File("");
			if(returnVal == JFileChooser.APPROVE_OPTION) {
				fileIn= chooser.getSelectedFile();
				currentPath = fileIn.getAbsolutePath();
			}
		}
		if ("setSaveFolder".equals(e.getActionCommand())){
			JFileChooser chooser = new JFileChooser(savePath);
			chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
			int returnVal = chooser.showOpenDialog(this);
			if(returnVal == JFileChooser.APPROVE_OPTION) {
				File selected = chooser.getSelectedFile();
				savePath = selected.getAbsolutePath();
			}
		}
		if ("setSaveImageFolder".equals(e.getActionCommand())){
			JFileChooser chooser = new JFileChooser(imageSavePath);
			chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
			int returnVal = chooser.showOpenDialog(this);
			if(returnVal == JFileChooser.APPROVE_OPTION) {
				File selected = chooser.getSelectedFile();
				imageSavePath = selected.getAbsolutePath();
			}
		}
		
		if ("analyze".equals(e.getActionCommand())) {
			//JOptionPane.showMessageDialog(this,"Nappia painettiin.");
			
			
			try{
				softThreshold = Double.parseDouble(softThresholdInput.getText());
				boneThreshold = Double.parseDouble(boneThresholdInput.getText());
				marrowThreshold = Double.parseDouble(marrowThresholdInput.getText());
				areaThreshold = Double.parseDouble(areaThresholdInput.getText());
				BMDThreshold = Double.parseDouble(BMDThresholdInput.getText());
				sectorWidth = 10;
				File[] listOfFiles;
				File tempFile = new File(currentPath);
				listOfFiles = tempFile.listFiles();
				details = new ImageAndAnalysisDetails(Double.parseDouble(tScaleCoeff.getText()),Double.parseDouble(tScaleConst.getText())
														,softThreshold, boneThreshold,marrowThreshold,areaThreshold,BMDThreshold
														,filterSize,softFilterSize, sectorWidth, imageSavePath
														, mRoiDet,dBoneSite,eOn,mOn,dOn, stOn, cOn,ukkOn,dicomOn,imOn);
				AnalysisThread analysisThread = new AnalysisThread(this,listOfFiles,showFigure, details);
				Thread anaThread = new Thread(analysisThread,"analysisThread");
				anaThread.start();	//All of the analysis needs to be done within this thread from hereafter
			}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		}
		
		/*
		if ("erode".equals(e.getActionCommand())) {
			//JOptionPane.showMessageDialog(this,"Nappia painettiin.");
			  try{
				  showFigure.drawImage(analyzeRoi.peeledROI,tiedot.PicMatrixX, tiedot.PicMatrixY,roi.minimum,roi.maximum);
				}catch (Exception err){System.err.println("Error: "+err.getMessage());}
			}
		if ("save".equals(e.getActionCommand())) {
			//JOptionPane.showMessageDialog(this,"Nappia painettiin.");
				
					
			}
			*/
	}
	
	public class PatientDetails{
		public String PatName;
		public String PatBirth;
		public String PatID;
		public String MeasDate;
		public String ObjLen;
		public String VoxelSize;
		public PatientDetails(ReadStratecFile tiedot){
			PatName = tiedot.PatName;
			PatBirth = Long.toString(tiedot.PatBirth);
			PatID = tiedot.PatID;
			MeasDate = Long.toString(tiedot.MeasDate);
			ObjLen = Double.toString(tiedot.ObjLen);
			VoxelSize = Double.toString(tiedot.VoxelSize);
		}
		public PatientDetails(DicomReader tiedot){
			PatName = new String(tiedot.patientName);
			PatBirth = new String(tiedot.patientBirthDate);
			PatID = new String(tiedot.patientID);
			MeasDate = new String(tiedot.acquisitionDate);
			ObjLen = new String(tiedot.imageComments);
			String pSpace =new String(tiedot.pixelSpacing);			
			VoxelSize = pSpace.substring(0,pSpace.indexOf(new String("\\"))-1);
		}
	}
	
	public class AnalysisThread implements Runnable{
		File[] listOfFiles;
		int sectorWidth;
		DrawImage showFigure;
		MidShaft_pQCT_analysis mainProgram;
		SelectROI roi;
		ImageAndAnalysisDetails details;
		ReadStratecFile readStratecFile;
		ScaledImageData scaledImageData;
		AnalyzeROI analyzeRoi;
		PatientDetails pDetails;
		DicomReader dicomReader;
		String imageSavePath;
		
		public AnalysisThread(MidShaft_pQCT_analysis mainProgramIn, File[] listOfFilesIn,DrawImage showFigureIn,
				ImageAndAnalysisDetails detailsIn){
				details = detailsIn;


			listOfFiles = listOfFilesIn;
//			sectorWidth = sectorWidthIn;
			showFigure = showFigureIn;
			mainProgram = mainProgramIn;

		}
		public void run(){
			
			for (int i = 0;i <listOfFiles.length; i++) {//i< 1;i++){//
				 if (listOfFiles[i].isFile()) {
					try{
							fileIn = listOfFiles[i];
						  BufferedInputStream tied = new BufferedInputStream(new FileInputStream(fileIn));
						  DataInputStream dataFile = new DataInputStream(tied);
							if (!dataFile.markSupported()) {
								throw new RuntimeException("Mark/Reset not supported!");
							}
						  //FilterInputStream tiedosto = new FilterInputStream(tied);
						  if (!dicomOn){
							  readStratecFile = new ReadStratecFile(dataFile);		//Read Data from file
							  scaledImageData = new ScaledImageData(readStratecFile, details);	//Scale and 3x3 median filter the data
							  pDetails = new PatientDetails(readStratecFile);
						  }else{
							dataFile.close();
							dicomReader = new DicomReader(fileIn);
							scaledImageData = new ScaledImageData(dicomReader, details);	//Scale and 3x3 median filter the data
							pDetails = new PatientDetails(dicomReader);
						  }
						  showFigure.drawImage(scaledImageData.scaledImage,scaledImageData.width, scaledImageData.height, scaledImageData.minimum, scaledImageData.maximum);
						   //Draw image on screen
						   //System.out.println("Figure on screen "+scaledImageData.minimum+" "+scaledImageData.maximum);
						  //USE ROI for selecting what to analyze --cortex, marrow (trabecular bone in the future)
						  roi = new SelectROI(scaledImageData, details,showFigure,fileIn.getName());
						  //new Thread(roi).start();'
						  //System.out.println("ROI ");
						  Thread roiThread = new Thread(roi,"roiThread");
						  roiThread.start();	//All of the analysis needs to be done within this thread from hereafter
						  roiThread.join();
						  //Begin analysis
						  //System.out.println("Info");
						  mainProgram.writeInfo(fileIn,pDetails,roi,details.sectorWidth);
						  //System.out.println("CoAna");
						  if (cOn){//cortical analysis
							 CorticalAnalysis cortAnalysis =new CorticalAnalysis(roi);
							 mainProgram.writeCortical(fileIn,cortAnalysis);
						  }
						  //System.out.println("MedAna");
						  if (mOn){//medullary analysis
							MedullaryAnalysis medAnalysis =new MedullaryAnalysis(roi);
							mainProgram.writeMedullary(fileIn,medAnalysis);
							}
						  //System.out.println("SoftAna");
						  if (stOn){//Soft tissue analysis
							SoftTissueAnalysis softTissueAnalysis = new SoftTissueAnalysis(roi);
							mainProgram.writeSoft(fileIn,softTissueAnalysis);
						  }
						  //System.out.println("DistAna");
						  if (dOn){//Distribution analysis
							  analyzeRoi = new AnalyzeROI(roi,details);
							  showFigure.drawImage(roi.cortexROI,scaledImageData.width, scaledImageData.height,roi.minimum,roi.maximum);
							  mainProgram.writeDistribution(fileIn, analyzeRoi,details.sectorWidth);
						  }
						  
                    }catch (Exception err){System.err.println("Error: "+err.getMessage());}
				  mainProgram.writeNewLine();
			  }
			}
		}
	}
	void writeHeader(int sectorWidth){	//Write header line when using a given result file for a first time
		try{//Write results to a file, use try - catch error handling
			FileWriter resultFile = new FileWriter(savePath+"/"+tNimi.getText(),true);	//Open for appending
			BufferedWriter writer = new BufferedWriter(resultFile);
			writer.write("Filename\t"
							+"Name\t"
							+"BirthDate\t"
							+"ID\t"
							+"Measurement date\t"
							+"Bone length [mm]\t"
							+"Pixel size [mm]\t"
							+"Leg [0 = right, 1 = left]\t");
							
							if (cOn){
			writer.write("CoD [mg/cm3]\t"
							+"CoA [mm2]\t"
							+"ToA [mm2]\t"
							+"BSI [mg2/dm4]\t"
							+"SSI [mm3]\t"
							+"dwIPo[mm4]\t"
							+"dwIMax [mm4]\t"
							+"dwIMin[mm4]\t"
							+"alfa [deg]\t");
							}
							if (mOn){
			writer.write("med vBMD [mg/cm3]\t"
							+"med Area [mm2]\t"
							+"kernel vBMD [mg/cm3]\t"
							+"kernel Area [mm2]\t"
							+"ring vBMD [mg/cm3]\t"
							+"ring Area [mm2]\t");
			for (int pp = 0;pp<10;pp++){
								writer.write("Concentric marrowRing "+pp+" [mg/cm$^3$]\t");
							}
							}
							if (stOn){
			writer.write(	"SoftTissueDensity [mg/cm3]\t"
							+"SoftTissueArea [mm2]\t"
							+"MuscleDensity [mg/cm3]\t"
							+"MuscleArea [mm2]\t"
							+"SubcutaneousFatDensity [mg/cm3]\t"
							+"SubcutaneousFatArea [mm2]\t"
							+"FatPercentage [%]\t"
							+"DensityWeightedFatPercentage [%]\t");
							}
			
				if (dOn){
					for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
						writer.write(pp*sectorWidth+" - "+((pp+1)*sectorWidth)+" endocortical radius [mm]\t");
					}
					for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
						writer.write(pp*sectorWidth+" - "+((pp+1)*sectorWidth)+" pericortical radius [mm]\t");
					}
				
					for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
						writer.write(pp*sectorWidth+" - "+((pp+1)*sectorWidth)+" endocortical radius after peeling [mm]\t");
					}
					for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
						writer.write(pp*sectorWidth+" - "+((pp+1)*sectorWidth)+" pericortical radius after peeling [mm]\t");
					}
					
					//Cortex BMD values			
					for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
						writer.write(pp*sectorWidth+" - "+((pp+1)*sectorWidth)+" endocortical vBMD [mg/cm3]\t");
					}
					for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
						writer.write(pp*sectorWidth+" - "+((pp+1)*sectorWidth)+" midcortical vBMD [mg/cm3]\t");
					}
					for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
						writer.write(pp*sectorWidth+" - "+((pp+1)*sectorWidth)+" pericortical vBMD [mg/cm3]\t");
					}
				}
			
			//Close the output stream
			writer.write("\n");
			writer.close();
		}catch (Exception err){//Catch exception if any
			System.err.println("Error: " + err.getMessage());
			marrowCenter.setText("Couldn't save");
		}
	}
	void writeNewLine(){
		try{
			FileWriter resultFile = new FileWriter(savePath+"/"+tNimi.getText(),true);	//Open for appending
			BufferedWriter writer = new BufferedWriter(resultFile);
			writer.write("\n");
			//Close the output stream
			writer.close();
		}catch (Exception err){//Catch exception if any
			System.err.println("Error: " + err.getMessage());
			marrowCenter.setText("Couldn't save");
		}
	}	

	void writeInfo(File fileIn,PatientDetails tiedot,SelectROI roi, int sectorWidth){
		try{
			File test = new File(savePath+"/"+tNimi.getText());
			if (test.exists() == false){writeHeader(sectorWidth);}
			FileWriter resultFile = new FileWriter(savePath+"/"+tNimi.getText(),true);	//Open for appending
			BufferedWriter writer = new BufferedWriter(resultFile);
			writer.write(fileIn.getName()+"\t"
							+tiedot.PatName+"\t"
							+tiedot.PatBirth+"\t"
							+tiedot.PatID+"\t"
							+tiedot.MeasDate+"\t"
							+tiedot.ObjLen+"\t"
							+tiedot.VoxelSize+"\t"
							+roi.leg+"\t");
			writer.close();//Close the output stream
		}catch (Exception err){//Catch exception if any
			System.err.println("Error: " + err.getMessage());
			marrowCenter.setText("Couldn't save");
		}
	}
	
	void writeMedullary(File fileIn,MedullaryAnalysis medAnalysis){
		try{
			FileWriter resultFile = new FileWriter(savePath+"/"+tNimi.getText(),true);	//Open for appending
			BufferedWriter writer = new BufferedWriter(resultFile);
			writer.write(medAnalysis.medBMD+"\t"
							+medAnalysis.medAREA+"\t"
							+medAnalysis.kernelBMD+"\t"
							+medAnalysis.kernelAREA+"\t"
							+medAnalysis.ringBMD+"\t"
							+medAnalysis.ringAREA+"\t");
			for (int i = 9;i>=0;i--){
				writer.write(medAnalysis.radialBMD[i]+"\t");
			}
			writer.close();	//Close the output stream
		}catch (Exception err){//Catch exception if any
			System.err.println("Error: " + err.getMessage());
			marrowCenter.setText("Couldn't save");
		}
	}
	void writeCortical(File fileIn,CorticalAnalysis cortAnalysis){
			try{
			FileWriter resultFile = new FileWriter(savePath+"/"+tNimi.getText(),true);	//Open for appending
			BufferedWriter writer = new BufferedWriter(resultFile);
			writer.write(cortAnalysis.BMD+"\t"
							+cortAnalysis.AREA+"\t"
							+cortAnalysis.ToA+"\t"
							+cortAnalysis.BMD*cortAnalysis.BMD*cortAnalysis.AREA/1000000.0+"\t"
							+cortAnalysis.SSI+"\t"
							+cortAnalysis.dwIPo+"\t"
							+cortAnalysis.dwIMax+"\t"
							+cortAnalysis.dwIMin+"\t"
							+cortAnalysis.alfa*180.0/Math.PI+"\t");
			writer.close();	//Close the output stream
		}catch (Exception err){//Catch exception if any
			System.err.println("Error: " + err.getMessage());
			marrowCenter.setText("Couldn't save");
		}
	}
	void writeSoft(File fileIn,SoftTissueAnalysis softTissueAnalysis){
			try{
			FileWriter resultFile = new FileWriter(savePath+"/"+tNimi.getText(),true);	//Open for appending
			BufferedWriter writer = new BufferedWriter(resultFile);
			writer.write(softTissueAnalysis.stBMD+"\t"
							+softTissueAnalysis.stAREA+"\t"
							+softTissueAnalysis.muscleBMD+"\t"
							+softTissueAnalysis.muscleAREA+"\t"
							+softTissueAnalysis.subBMD+"\t"
							+softTissueAnalysis.subAREA+"\t"
							+softTissueAnalysis.fatPercentage+"\t"
							+softTissueAnalysis.fatPercentageWeighted+"\t");
			writer.close();	//Close the output stream
		}catch (Exception err){//Catch exception if any
			System.err.println("Error: " + err.getMessage());
			marrowCenter.setText("Couldn't save");
		}
	}
	void writeDistribution(File fileIn,AnalyzeROI analyzeRoi, int sectorWidth){
		try{
			FileWriter resultFile = new FileWriter(savePath+"/"+tNimi.getText(),true);	//Open for appending
			BufferedWriter writer = new BufferedWriter(resultFile);
			
			for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
				writer.write(analyzeRoi.endocorticalRadii[pp]+"\t");
			}
			for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
				writer.write(analyzeRoi.pericorticalRadii[pp]+"\t");
			}
		
			for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
				writer.write(analyzeRoi.peeledEndocorticalRadii[pp]+"\t");
			}
			for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
				writer.write(analyzeRoi.peeledPericorticalRadii[pp]+"\t");
			}
			
			//Cortex BMD values			
			for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
				writer.write(analyzeRoi.endoCorticalBMDs[pp]+"\t");
			}
			for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
				writer.write(analyzeRoi.midCorticalBMDs[pp]+"\t");
			}
			for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
				writer.write(analyzeRoi.periCorticalBMDs[pp]+"\t");
			}
			//Close the output stream
			writer.close();
		}catch (Exception err){//Catch exception if any
			System.err.println("Error: " + err.getMessage());
			marrowCenter.setText("Couldn't save");
		}
	}
	public static void main(String[] args){
		javax.swing.SwingUtilities.invokeLater(new Runnable() {
			public void run(){
				doAndShowGUI();
			}
		});
	
		

	}
	    /* Listens to the check boxes. */
    public void itemStateChanged(ItemEvent e) {
        Object source = e.getItemSelectable();
		if (source == manualRoiDetermination) {
			if (e.getStateChange() == ItemEvent.DESELECTED) {
				mRoiDet = false;
			}
			if (e.getStateChange() == ItemEvent.SELECTED) {
				mRoiDet = true;
			}
			System.out.println("ManualROI "+mRoiDet);
		}
		if (source == erodeOn) {
			if (e.getStateChange() == ItemEvent.DESELECTED) {
				eOn = false;
			}
			if (e.getStateChange() == ItemEvent.SELECTED) {
				eOn = true;
			}
			System.out.println("Erode "+eOn);
		}
		if (source == marrowAna) {
			if (e.getStateChange() == ItemEvent.DESELECTED) {
				mOn = false;
			}
			if (e.getStateChange() == ItemEvent.SELECTED) {
				mOn = true;
			}
			System.out.println("MarrowAnalysis "+mOn);
		}
		if (source == distributionAna) {
			if (e.getStateChange() == ItemEvent.DESELECTED) {
				dOn = false;
			}
			if (e.getStateChange() == ItemEvent.SELECTED) {
				dOn = true;
			}
			System.out.println("DistributionAnalysis "+dOn);
		}
		if (source == softAna) {
			if (e.getStateChange() == ItemEvent.DESELECTED) {
				stOn = false;
			}
			if (e.getStateChange() == ItemEvent.SELECTED) {
				stOn = true;
			}
			System.out.println("Soft Tissue Analysis "+stOn);
		}
		if (source == cortAna) {
			if (e.getStateChange() == ItemEvent.DESELECTED) {
				cOn = false;
			}
			if (e.getStateChange() == ItemEvent.SELECTED) {
				cOn = true;
			}
			System.out.println("Cortical Analysis "+cOn);
		}
		
		if (source == ukkSpecial) {
			if (e.getStateChange() == ItemEvent.DESELECTED) {
				ukkOn = false;
			}
			if (e.getStateChange() == ItemEvent.SELECTED) {
				ukkOn = true;
			}
			System.out.println("UKK special "+ukkOn);
		}
		if (source == distalBone) {
			if (e.getStateChange() == ItemEvent.DESELECTED) {
				dBoneSite = false;
				boneThreshold  = 550;
				boneThresholdInput.setText(oneDec.format(boneThreshold));

				
			}
			if (e.getStateChange() == ItemEvent.SELECTED) {
				dBoneSite = true;
				boneThreshold  = 169;
				boneThresholdInput.setText(oneDec.format(boneThreshold));
			}
			System.out.println("ManualROI "+mRoiDet);
		}
		if (source == dicomAna) {
			if (e.getStateChange() == ItemEvent.DESELECTED) {
				dicomOn = false;
			}
			if (e.getStateChange() == ItemEvent.SELECTED) {
				dicomOn = true;
			}
			System.out.println("Dicom files analyzed "+dicomOn);
		}
		
		if (source == imSave) {
			if (e.getStateChange() == ItemEvent.DESELECTED) {
				imOn = false;
			}
			if (e.getStateChange() == ItemEvent.SELECTED) {
				imOn = true;
			}
			System.out.println("Visual inspection images saved "+imOn);
		}
    }
}

	
