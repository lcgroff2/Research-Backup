/**
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */

package biotransformer.utils;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;
import org.codehaus.jackson.JsonParseException;
import org.codehaus.jackson.map.JsonMappingException;
import org.json.simple.parser.ParseException;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.btransformers.Biotransformer;
import biotransformer.btransformers.Cyp450BTransformer;
import biotransformer.btransformers.ECBasedBTransformer;
import biotransformer.btransformers.HGutBTransformer;
import biotransformer.btransformers.Phase2BTransformer;
import biotransformer.transformation.Biotransformation;
import biotransformer.transformation.MetabolicReaction;
import biotransformer.transformation.MetabolicPathway.MPathwayName;
import biotransformer.utils.ChemicalClassFinder.ChemicalClassName;
import predicition.P2Filter;

public class HumanSuperBioTransformer {
	
	protected ECBasedBTransformer ecb		 		= new ECBasedBTransformer(BioSystemName.HUMAN);
	protected Cyp450BTransformer cyb 				= new Cyp450BTransformer(BioSystemName.HUMAN);
	protected HGutBTransformer hgb 					= new HGutBTransformer();
	protected Phase2BTransformer p2b 				= new Phase2BTransformer(BioSystemName.HUMAN);
	protected LinkedHashMap<String, LinkedHashMap<String, String>> compoundDictionary		
													= new LinkedHashMap<String, LinkedHashMap<String, String>>();
	
	protected P2Filter p2filter 					= new P2Filter();
	protected LinkedHashMap<String, MetabolicReaction> combinedReactionsHash							
													= new LinkedHashMap<String, MetabolicReaction>();
	
	public SmilesParser smiParser					= ecb.getSmiParser();
	public SmilesGenerator smiGen 		= new SmilesGenerator().isomeric();
	
	public HumanSuperBioTransformer() throws IOException, ParseException, CDKException {

		for(Map.Entry<String, MetabolicReaction> m : this.ecb.reactionsHash.entrySet()){
			if(! this.combinedReactionsHash.containsKey(m.getKey())){
				this.combinedReactionsHash.put(m.getKey(), m.getValue());
			}
		}
		
		for(Map.Entry<String, MetabolicReaction> n : this.cyb.reactionsHash.entrySet()){
			if(! this.combinedReactionsHash.containsKey(n.getKey())){
				this.combinedReactionsHash.put(n.getKey(), n.getValue());
			}
		}
//		System.err.println("No. of reactions: " + this.combinedReactionsHash.size());
//		System.err.println(this.combinedReactionsHash.containsKey("STEROL_16_HYDROXYLATION_PATTERN2"));
		
		for(Map.Entry<String, MetabolicReaction> p : this.hgb.reactionsHash.entrySet()){
			if(! this.combinedReactionsHash.containsKey(p.getKey())){
				this.combinedReactionsHash.put(p.getKey(), p.getValue());
			}
		}
//		System.err.println("No. of reactions: " + this.combinedReactionsHash.size());
		
		for(Map.Entry<String, MetabolicReaction> o : this.p2b.reactionsHash.entrySet()){
			if(! this.combinedReactionsHash.containsKey(o.getKey())){
				this.combinedReactionsHash.put(o.getKey(), o.getValue());
			}
		}
//		System.err.println("No. of reactions: " + this.combinedReactionsHash.size());
		
//		System.err.println("Combined reactionsHash");
//		System.err.println(this.combinedReactionsHash);
//		System.err.println(this.combinedReactionsHash.size());
//		System.out.println(this.cyb.reactionsHash);
		
//		System.err.println("No. of reactions: " + this.combinedReactionsHash.size());
//		for(Map.Entry<String, MetabolicReaction> m : this.combinedReactionsHash.entrySet()){
//			System.out.println(m.getKey());
//		}
//		System.err.println(this.combinedReactionsHash.containsKey("EAWAG_RULE_BT0008"));
	}

	
	public ArrayList<Biotransformation> simulateHumanSuperbioMetabolism(IAtomContainer target) throws Exception{
		 ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		 
//		 	System.out.println("Is Biotransformer valid: " + ChemStructureExplorer.isBioTransformerValid(target));
			if(ChemStructureExplorer.isBioTransformerValid(target)){
				 IAtomContainer molecule = ChemStructureManipulator.standardizeMoleculeWithCopy(target);
//				 AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
				 IAtomContainerSet products = DefaultChemObjectBuilder
							.getInstance().newInstance(IAtomContainerSet.class);
				 
				 IAtomContainerSet phaseIISubstrates = DefaultChemObjectBuilder
							.getInstance().newInstance(IAtomContainerSet.class);
				 
				 IAtomContainerSet phaseIINonSubstrates = DefaultChemObjectBuilder
							.getInstance().newInstance(IAtomContainerSet.class);
				 
				 products.addAtomContainer(molecule);
				 
				 ChemicalClassName chemClassName = ChemicalClassFinder.findChemicalClass(molecule);
				 
				 
//				if(!(ChemicalClassFinder.isEtherLipid(molecule) || ChemicalClassFinder.isGlyceroLipid(molecule) || ChemicalClassFinder.isGlycerophosphoLipid(molecule) ||
//						ChemicalClassFinder.isSphingoLipid(molecule))) {
					if(!(chemClassName == ChemicalClassName.ETHER_LIPID || chemClassName == ChemicalClassName.GLYCEROLIPID || chemClassName == ChemicalClassName.GLYCEROPHOSPHOLIPID ||
							chemClassName == ChemicalClassName.SPHINGOLIPID || chemClassName == ChemicalClassName.GLYCEROL_3_PHOSPHATE_INOSITOL )) {	
//					 System.out.println("Predicting EC  metabolism round 1");
					/*
					 *  Apply ECBased
					 */		 
					
					System.out.println("\n\n===========================================");
					System.out.println("Predicting EC metabolism");
					System.out.println("===========================================\n\n");
					System.out.println("Smiles before EC-based simulation: " + this.smiGen.create(molecule));
					 biotransformations.addAll(this.ecb.simulateECBasedPhaseIMetabolismChain(molecule, true, true, 1, 0.1));
					 products.add(this.ecb.extractAtomContainer(biotransformations));
					 products = ChemStructureExplorer.uniquefy(products);
					 System.out.println("Number of EC-based biotransformations after first pass: " + biotransformations.size());
//					 System.out.println("Predicting CYP450 metabolism");
					/*
					 *  Apply CYP450 metabolism
					 */	 
					 
					System.out.println("\n\n===========================================");
					System.out.println("Predicting CYP450 metabolism");
					System.out.println("===========================================\n\n");
					 ArrayList<Biotransformation> cyp450Biots = new  ArrayList<Biotransformation>();
					 for(IAtomContainer met : products.atomContainers()){
//						 System.out.println("predicting CYP450 metabolites for: " + this.ecb.smiGen.create(met));
						 cyp450Biots.addAll(this.cyb.predictCyp450Biotransformations(met, true, true, 0.1));
					 }
					 IAtomContainerSet cyp450Prods = this.ecb.extractAtomContainer(cyp450Biots);
					 products.add(cyp450Prods);
					 products = ChemStructureExplorer.uniquefy(products);
					 biotransformations.addAll(cyp450Biots);
//					 System.out.println("Number of CYP450 biotransformations: " + cyp450Biots.size());
					 System.out.println("Number of CYP450 products: " + cyp450Prods.getAtomContainerCount());
//					 for(IAtomContainer a : cyp450Prods.atomContainers()){
//						 System.out.println(this.smiGen.create(a));
//					 }
					 
					 System.out.println("products: " + products.getAtomContainerCount());
					 
					 
					/*
					 *  Apply ECBased
					 *  Some products of CYP450 might be unstable, such as epoxides, and will be further transformed.
					 *  
					 *  !!!!!!!! MAKE SURE THESE DO NOT INCLUDE HYDROLYSIS REATIONS AND SOME OTHERS. THAT OFTEN DO NOT OCCUR AFTER CYP450
					 */		 
					 
					 LinkedHashMap<String, IAtomContainerSet> partitionedMolecules = this.p2filter.partitionSetForPhaseIIMetabolism(products);
					 
					 phaseIISubstrates.add(partitionedMolecules.get("phaseIISubstrates"));
					 phaseIINonSubstrates.add(partitionedMolecules.get("phaseIINonSubstrates"));
					 
					 System.out.println("phaseIISubstrates: " + phaseIISubstrates.getAtomContainerCount());
					 
					 if	(phaseIINonSubstrates.getAtomContainerCount()>0){
						System.out.println("\n\n===========================================");
						System.out.println("Predicting EC metabolism - 2nd pass");
						System.out.println("===========================================\n\n");
						 ArrayList<Biotransformation> ecBiotsSecondPass = new  ArrayList<Biotransformation>();
						 ecBiotsSecondPass = this.ecb.simulateECBasedPhaseIMetabolismChain(partitionedMolecules.get("phaseIINonSubstrates"), true, true, 1, 0.1);
						 biotransformations.addAll(ecBiotsSecondPass);
						 IAtomContainerSet ecBiotsSecondPassProducts = this.ecb.extractAtomContainer(ecBiotsSecondPass);					 
						 products.add(ecBiotsSecondPassProducts);
						 products = ChemStructureExplorer.uniquefy(products);
						 
						 /*
						  *   Only the ones suitable for phaseII will land into the gut
						  */
						 
						 LinkedHashMap<String, IAtomContainerSet> partitionedMoleculesAfterEC2 = this.p2filter.partitionSetForPhaseIIMetabolism(ecBiotsSecondPassProducts);
//						 System.out.print("partitionedMoleculesAfterEC2");
//						 System.out.print("phaseIISubstrates" + partitionedMoleculesAfterEC2.get("phaseIISubstrates").getAtomContainerCount());
//						 System.out.print("phaseIINonSubstrates" + partitionedMoleculesAfterEC2.get("phaseIINonSubstrates").getAtomContainerCount());
						 
						 phaseIISubstrates.add(partitionedMoleculesAfterEC2.get("phaseIISubstrates"));
						 
						 phaseIINonSubstrates.add(partitionedMoleculesAfterEC2.get("phaseIINonSubstrates"));	
//						 System.out.println("Number of EC-based biotransformations during second pass: " + ecBiotsSecondPass.size());
//						 System.out.println("Number of EC-based metabolites during second pass: " + ecBiotsSecondPassProducts.getAtomContainerCount());
					 }
					 

					 System.out.println("Predicting human gut metabolism of " + products.getAtomContainerCount() + " metabolites");
					/*
					 *  Apply Human gut metabolism
					 */
					 
//					 for(IAtomContainer atc : products.atomContainers()){
//						 System.out.println(this.ecb.smiGen.create(atc));
//					 }
					 
					 IAtomContainerSet hGutSubstrates = phaseIISubstrates;
					 
					 for(IAtomContainer a : phaseIINonSubstrates.atomContainers()){
						 /*
						  * add some large molecules, such as tannins, which can be degraded by bacteria
						  */
						 if(ChemStructureExplorer.getMajorIsotopeMass(a) >= 900.0){
							 hGutSubstrates.addAtomContainer(a);
						 }
					 }
					 
//					 System.out.println(hGutSubstrates == null);
//					 System.out.println(hGutSubstrates.getAtomContainerCount());
					 ArrayList<Biotransformation> hGutBiots = new  ArrayList<Biotransformation>();
					 hGutBiots.addAll(this.hgb.simulateGutMicrobialMetabolismHydrolysisAndReduction(hGutSubstrates, true, true, 8, 0.1));
					 System.out.println("Number of human gut biotransformations: " + hGutBiots.size());
					 biotransformations.addAll(hGutBiots);
					 IAtomContainerSet hGutProducts = this.ecb.extractAtomContainer(hGutBiots);
					 System.out.println("Number of human gut metabolites: " + hGutProducts.getAtomContainerCount());
					 products.add(hGutProducts);
//					 products = ChemStructureExplorer.uniquefy(products);
					 
					 LinkedHashMap<String, IAtomContainerSet> partitionedMoleculesAfterHGut = this.p2filter.partitionSetForPhaseIIMetabolism(hGutProducts);
					 
					 phaseIISubstrates.add(partitionedMoleculesAfterHGut.get("phaseIISubstrates"));
					 phaseIINonSubstrates.add(partitionedMoleculesAfterHGut.get("phaseIINonSubstrates"));
					 
					 phaseIISubstrates = ChemStructureExplorer.uniquefy(phaseIISubstrates);
					 
//					 System.out.println("Predicting phase II metabolism of " + phaseIISubstrates.getAtomContainerCount() + " out of " + products.getAtomContainerCount());
					 
//					 for(IAtomContainer atomc : phaseIISubstrates.atomContainers()){
//						 System.out.println(this.smiGen.create(atomc));
//					 }
					 
					 
					 if(phaseIISubstrates.getAtomContainerCount()>0){
						 
//						 IAtomContainerSet phase2Cadidates = this.p2filter.returnFilteredPhaseIICandidates(products);
						 
//						 System.out.println("Predicting phase II metabolism of " + phaseIISubstrates.getAtomContainerCount() + " out of " + products.getAtomContainerCount());
						
//						 for(IAtomContainer atc : phaseIISubstrates.atomContainers()){
//							 System.out.println(this.ecb.smiGen.create(atc));
//						 }
						 
						 /*
						 *  Apply Phase II metabolism
						 */	 
						System.out.println("Predicting Phase 2 metabolism for " + phaseIISubstrates.getAtomContainerCount());
						System.out.println("\n\n===========================================");
						System.out.println("Predicting human phase 2 metabolism");
						System.out.println("===========================================\n\n");
						 ArrayList<Biotransformation> phaseIIBiots = new  ArrayList<Biotransformation>();
						 phaseIIBiots.addAll(this.p2b.applyPhase2TransformationsChainAndReturnBiotransformations(phaseIISubstrates, true, true, true, 1, 0.1));
						 products.add(this.ecb.extractAtomContainer(phaseIIBiots));
						 products = ChemStructureExplorer.uniquefy(products);
						 biotransformations.addAll(phaseIIBiots);
						 System.out.println("Number of PhaseII biotransformations: " + phaseIIBiots.size() + "\n\n\n");				
					 }

				}
					
				else {
//					System.out.println("\n\n===========================================");
//					System.out.println("Predicting EC metabolism");
//					System.out.println("===========================================\n\n");
//					ArrayList<Biotransformation> ecBiots = this.ecb.simulateECBasedMetabolismChain(molecule, true, true, 1, 0.1);
//					biotransformations.addAll(ecBiots);
//					System.out.println("Number of EC biotransformations: " + ecBiots.size() + "\n\n\n");

					
//						 System.out.println("Predicting EC  metabolism round 1");
					/*
					 *  Apply ECBased
					 */		 
					
					System.out.println("\n\n===========================================");
					System.out.println("Predicting EC metabolism");
					System.out.println("===========================================\n\n");
					System.out.println("Smiles before EC-based simulation: " + this.smiGen.create(molecule));
					 biotransformations.addAll(this.ecb.simulateECBasedPhaseIMetabolismChain(molecule, true, true, 1, 0.1));
					 products.add(this.ecb.extractAtomContainer(biotransformations));
					 products = ChemStructureExplorer.uniquefy(products);
					 System.out.println("Number of EC-based biotransformations after first pass: " + biotransformations.size());

					 /*
					 *  Apply ECBased
					 *  Some products of CYP450 might be unstable, such as epoxides, and will be further transformed.
					 *  
					 *  !!!!!!!! MAKE SURE THESE DO NOT INCLUDE HYDROLYSIS REATIONS AND SOME OTHERS. THAT OFTEN DO NOT OCCUR AFTER CYP450
					 */		 
					 
					 LinkedHashMap<String, IAtomContainerSet> partitionedMolecules = this.p2filter.partitionSetForPhaseIIMetabolism(products);
					 
					 phaseIISubstrates.add(partitionedMolecules.get("phaseIISubstrates"));
					 phaseIINonSubstrates.add(partitionedMolecules.get("phaseIINonSubstrates"));
					 
					 System.out.println("phaseIISubstrates: " + phaseIISubstrates.getAtomContainerCount());
					 
//					 if	(phaseIINonSubstrates.getAtomContainerCount()>0){
//						System.out.println("\n\n===========================================");
//						System.out.println("Predicting EC metabolism - 2nd pass");
//						System.out.println("===========================================\n\n");
//						 ArrayList<Biotransformation> ecBiotsSecondPass = new  ArrayList<Biotransformation>();
//						 ecBiotsSecondPass = this.ecb.simulateECBasedPhaseIMetabolismChain(partitionedMolecules.get("phaseIINonSubstrates"), true, true, 1, 0.1);
//						 biotransformations.addAll(ecBiotsSecondPass);
//						 IAtomContainerSet ecBiotsSecondPassProducts = this.ecb.extractAtomContainer(ecBiotsSecondPass);					 
//						 products.add(ecBiotsSecondPassProducts);
//						 products = ChemStructureExplorer.uniquefy(products);
//						 
//						 /*
//						  *   Only the ones suitable for phaseII will land into the gut
//						  */
//						 
//						 LinkedHashMap<String, IAtomContainerSet> partitionedMoleculesAfterEC2 = this.p2filter.partitionSetForPhaseIIMetabolism(ecBiotsSecondPassProducts);
////							 System.out.print("partitionedMoleculesAfterEC2");
////							 System.out.print("phaseIISubstrates" + partitionedMoleculesAfterEC2.get("phaseIISubstrates").getAtomContainerCount());
////							 System.out.print("phaseIINonSubstrates" + partitionedMoleculesAfterEC2.get("phaseIINonSubstrates").getAtomContainerCount());
//						 
//						 phaseIISubstrates.add(partitionedMoleculesAfterEC2.get("phaseIISubstrates"));
//						 
//						 phaseIINonSubstrates.add(partitionedMoleculesAfterEC2.get("phaseIINonSubstrates"));	
////							 System.out.println("Number of EC-based biotransformations during second pass: " + ecBiotsSecondPass.size());
////							 System.out.println("Number of EC-based metabolites during second pass: " + ecBiotsSecondPassProducts.getAtomContainerCount());
//					 }
					 

					 System.out.println("Predicting human gut metabolism of " + products.getAtomContainerCount() + " metabolites");
					/*
					 *  Apply Human gut metabolism
					 */
					 
//						 for(IAtomContainer atc : products.atomContainers()){
//							 System.out.println(this.ecb.smiGen.create(atc));
//						 }
					 
					 IAtomContainerSet hGutSubstrates = phaseIISubstrates;
					 
					 for(IAtomContainer a : phaseIINonSubstrates.atomContainers()){
						 /*
						  * add some large molecules, such as tannins, which can be degraded by bacteria
						  */
						 if(ChemStructureExplorer.getMajorIsotopeMass(a) >= 900.0){
							 hGutSubstrates.addAtomContainer(a);
						 }
					 }
					 
//					 System.out.println(hGutSubstrates == null);
//					 System.out.println(hGutSubstrates.getAtomContainerCount());
					 ArrayList<Biotransformation> hGutBiots = new  ArrayList<Biotransformation>();
					 hGutBiots.addAll(this.hgb.simulateGutMicrobialMetabolismHydrolysisAndReduction(hGutSubstrates, true, true, 1, 0.1));
					 System.out.println("Number of human gut biotransformations: " + hGutBiots.size());
					 biotransformations.addAll(hGutBiots);
					 IAtomContainerSet hGutProducts = this.ecb.extractAtomContainer(hGutBiots);
					 System.out.println("Number of human gut metabolites: " + hGutProducts.getAtomContainerCount());
					 products.add(hGutProducts);
//						 products = ChemStructureExplorer.uniquefy(products);
					 
					 LinkedHashMap<String, IAtomContainerSet> partitionedMoleculesAfterHGut = this.p2filter.partitionSetForPhaseIIMetabolism(hGutProducts);
					 
					 phaseIISubstrates.add(partitionedMoleculesAfterHGut.get("phaseIISubstrates"));
					 phaseIINonSubstrates.add(partitionedMoleculesAfterHGut.get("phaseIINonSubstrates"));
					 
					 phaseIISubstrates = ChemStructureExplorer.uniquefy(phaseIISubstrates);
					 
//						 System.out.println("Predicting phase II metabolism of " + phaseIISubstrates.getAtomContainerCount() + " out of " + products.getAtomContainerCount());
					 
//						 for(IAtomContainer atomc : phaseIISubstrates.atomContainers()){
//							 System.out.println(this.smiGen.create(atomc));
//						 }
					 
					 
					 if(phaseIISubstrates.getAtomContainerCount()>0){
						 
//							 IAtomContainerSet phase2Cadidates = this.p2filter.returnFilteredPhaseIICandidates(products);
						 
//							 System.out.println("Predicting phase II metabolism of " + phaseIISubstrates.getAtomContainerCount() + " out of " + products.getAtomContainerCount());
						
//							 for(IAtomContainer atc : phaseIISubstrates.atomContainers()){
//								 System.out.println(this.ecb.smiGen.create(atc));
//							 }
						 
						 /*
						 *  Apply Phase II metabolism
						 */	 
						System.out.println("Predicting Phase 2 metabolism for " + phaseIISubstrates.getAtomContainerCount() + " metabolites.");
						System.out.println("\n\n===========================================");
						System.out.println("Predicting human phase 2 metabolism");
						System.out.println("===========================================\n\n");
						 ArrayList<Biotransformation> phaseIIBiots = new  ArrayList<Biotransformation>();
						 phaseIIBiots.addAll(this.p2b.applyPhase2TransformationsChainAndReturnBiotransformations(phaseIISubstrates, true, true, true, 1, 0.1));
						 products.add(this.ecb.extractAtomContainer(phaseIIBiots));
						 products = ChemStructureExplorer.uniquefy(products);
						 biotransformations.addAll(phaseIIBiots);
						 System.out.println("Number of PhaseII biotransformations: " + Utilities.selectUniqueBiotransformations(phaseIIBiots).size() + "\n\n\n");				
					 }
				}

			}
		 
		IAtomContainerSet a = this.extractAtomContainer(biotransformations);
		System.out.println("Number of predicted compounds: " +  ChemStructureExplorer.uniquefy(a).getAtomContainerCount());
		 
		 System.out.println("Number of predicted biotransformations: " + biotransformations.size());
		 ArrayList<Biotransformation> uniqueBiotransformations = Utilities.selectUniqueBiotransformations(biotransformations);
		 System.out.println("Number of unique predicted biotransformations: " + uniqueBiotransformations.size());
		 return uniqueBiotransformations;
	}

	public void simulateHumanSuperbioMetabolismAndSaveToSDF(IAtomContainer molecule, String outputFileName, boolean annotate) throws Exception{
		ArrayList<Biotransformation> biotransformations = this.simulateHumanSuperbioMetabolism(molecule);

		this.ecb.saveBioTransformationProductsToSdf(Utilities.selectUniqueBiotransformations(biotransformations), outputFileName, this.combinedReactionsHash, annotate);
	}
		
	public void simulateHumanSuperbioMetabolismFromSDF(String inputFileName, boolean annotate) throws Exception {
		simulateHumanSuperbioMetabolismFromSDF(inputFileName, "data", annotate);
	}
	
	public void simulateHumanSuperbioMetabolismFromSDF(String inputFileName, String outputFolder, boolean annotate) throws Exception {
		
		IAtomContainerSet containers = FileUtilities.parseSdf(inputFileName);
		int nr = 0;
		if(!containers.isEmpty()){
			for(IAtomContainer molecule : containers.atomContainers()){
				try{
					String identifier = molecule.getProperty(CDKConstants.TITLE);
					if(identifier == null){
						identifier = molecule.getProperty("Name");
						if(identifier == null){
							identifier = molecule.getProperty("$MolName"); 
							if(identifier == null){
								identifier = molecule.getProperty("InChiKey");
								if(identifier == null){
									identifier = this.ecb.inchiGenFactory.getInChIGenerator(molecule).getInchiKey();
								}
							}

						}
					}
					identifier = identifier.replace(":", "-").replace("/", "_");
					this.simulateHumanSuperbioMetabolismAndSaveToSDF(molecule, outputFolder + "/" + identifier + "_BioT_sim_metabolites.sdf", annotate);

				}
				catch(Exception e){
					System.err.println("Could not predicted metabolism for molecule nr. " + nr);
					System.err.println(e.getLocalizedMessage());
				}
			}			
		}
		
	}
	
	public void  simulateHumanSuperbioMetabolismAndSaveToSDF(IAtomContainerSet containers, String outputFolder, boolean annotate) throws Exception {
		int nr = 0;
		if(!containers.isEmpty()){
			for(IAtomContainer molecule : containers.atomContainers()){

				try{
					String identifier = molecule.getProperty(CDKConstants.TITLE);
					if(identifier == null){
						identifier = molecule.getProperty("Name");
						if(identifier == null){
							identifier = molecule.getProperty("$MolName"); 
							if(identifier == null){
								identifier = molecule.getProperty("InChiKey");
								if(identifier == null){
									identifier = this.ecb.inchiGenFactory.getInChIGenerator(molecule).getInchiKey();
								}
							}

						}
					}
					identifier = identifier.replace(":", "-").replace("/", "_");
					this.simulateHumanSuperbioMetabolismAndSaveToSDF(molecule, outputFolder + "/" + identifier + "_BioT_sim_metabolites.sdf", annotate);

				}
				catch(Exception e){
					System.err.println("Could not predicted metabolism for molecule nr. " + nr);
					System.err.println(e.getLocalizedMessage());
				}

			}
		}		
	}
	
	
	
	
	public void simulateHumanSuperbioMetabolismAndSaveToCSV(IAtomContainer molecule, String outputFileName, boolean annotate) throws Exception{
		ArrayList<Biotransformation> biotransformations = this.simulateHumanSuperbioMetabolism(molecule);

		this.ecb.saveBioTransformationProductsToCSV(Utilities.selectUniqueBiotransformations(biotransformations), outputFileName, this.combinedReactionsHash, annotate);
	}
		
	
//	public void  simulateHumanMetabolismAndSaveToCSV(IAtomContainerSet containers, String outputFolder, boolean annotate) throws Exception {
//		int nr = 0;
//		if(!containers.isEmpty()){
//			for(IAtomContainer molecule : containers.atomContainers()){
//
//				try{
//					String identifier = molecule.getProperty(CDKConstants.TITLE);
//					if(identifier == null){
//						identifier = molecule.getProperty("Name");
//						if(identifier == null){
//							identifier = molecule.getProperty("$MolName"); 
//							if(identifier == null){
//								identifier = molecule.getProperty("InChiKey");
//								if(identifier == null){
//									identifier = this.ecb.inchiGenFactory.getInChIGenerator(molecule).getInchiKey();
//								}
//							}
//
//						}
//					}
//					identifier = identifier.replace(":", "-").replace("/", "_");
//					this.simulateHumanAndGutMicrobialMetabolismAndSaveToCSV(molecule, outputFolder + "/" + identifier + "_BioT_sim_metabolites.sdf", annotate);
//
//				}
//				catch(Exception e){
//					System.err.println("Could not predicted metabolism for molecule nr. " + nr);
//					System.err.println(e.getLocalizedMessage());
//				}
//
//			}
//		}		
//	}
	
	
	
	public ArrayList<Biotransformation>  simulateAllHumanMetabolism(IAtomContainerSet containers, double scoreThreshold) throws Exception {
		int nr = 0;
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		if(!containers.isEmpty()){
			for(IAtomContainer molecule : containers.atomContainers()){

				try{
					biotransformations.addAll(simulateOneStepAllHuman(molecule, scoreThreshold));
				}
				catch(Exception e){
					System.err.println("Could not predicted metabolism for molecule nr. " + nr);
					System.err.println(e.getLocalizedMessage());
				}

			}
		}
		return biotransformations;
	}
	

	
	
	public void simulateAllHumanMetabolismAndSavetoCSV(IAtomContainer molecule, String outputFileName, double scoreThreshold, boolean annotate) throws Exception {
		try{
			ArrayList<Biotransformation> biotransformations = simulateOneStepAllHuman(molecule, scoreThreshold);
			this.ecb.saveBioTransformationProductsToCSV(Utilities.selectUniqueBiotransformations(biotransformations), outputFileName, this.combinedReactionsHash, annotate);
		}
		catch(Exception e){
			System.err.println(e.getLocalizedMessage());
		}		
	}

	
	public void simulateAllHumanMetabolismAndSavetoSDF(IAtomContainer molecule, String outputFileName, double scoreThreshold, boolean annotate) throws Exception {
		try{
			ArrayList<Biotransformation> biotransformations = simulateOneStepAllHuman(molecule, scoreThreshold);
			this.ecb.saveBioTransformationProductsToSdf(Utilities.selectUniqueBiotransformations(biotransformations), outputFileName, this.combinedReactionsHash, annotate);
		}
		catch(Exception e){
			System.err.println(e.getLocalizedMessage());
		}		
	}
	
	
	
	
	public void simulateAllHumanMetabolismAndSavetoCSV(IAtomContainerSet containers, String outputFileName, double scoreThreshold, boolean annotate) throws Exception {
		try{
			ArrayList<Biotransformation> biotransformations = simulateAllHumanMetabolism(containers, scoreThreshold);
			this.ecb.saveBioTransformationProductsToCSV(Utilities.selectUniqueBiotransformations(biotransformations), outputFileName, this.combinedReactionsHash, annotate);
		}
		catch(Exception e){
			System.err.println(e.getLocalizedMessage());
		}		
	}

	
	public void simulateAllHumanMetabolismAndSavetoSDF(IAtomContainerSet containers, String outputFileName, double scoreThreshold, boolean annotate) throws Exception {
		try{
			ArrayList<Biotransformation> biotransformations = simulateAllHumanMetabolism(containers, scoreThreshold);
			this.ecb.saveBioTransformationProductsToSdf(Utilities.selectUniqueBiotransformations(biotransformations), outputFileName, this.combinedReactionsHash, annotate);
		}
		catch(Exception e){
			System.err.println(e.getLocalizedMessage());
		}		
	}
	

	public LinkedHashMap<String, IAtomContainerSet> partitionSetForPhaseIIMetabolism(IAtomContainerSet products) throws CDKException, Exception{
		return this.p2filter.partitionSetForPhaseIIMetabolism(products);
	}
	
	public void simulateHumanSuperbioMetabolismFromSDFtoSingleSDF(String inputFileName, String outputFileName, boolean annotate) throws Exception {
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
				
		int nr = 0;
		IAtomContainerSet containers = FileUtilities.parseSdf(inputFileName);
		if(!containers.isEmpty()){
			for(IAtomContainer molecule : containers.atomContainers()){
				nr++;
				System.out.println("molecule nr. " +nr);
				
				try{
					String identifier = molecule.getProperty(CDKConstants.TITLE);
					if(identifier == null){
						identifier = molecule.getProperty("Name");
						if(identifier == null){
							identifier = molecule.getProperty("$MolName"); 
							if(identifier == null){
								identifier = molecule.getProperty("InChiKey");
								if(identifier == null){
									identifier = this.ecb.inchiGenFactory.getInChIGenerator(molecule).getInchiKey();
								}
							}

						}
					}
					
					ArrayList<Biotransformation> bts = this.simulateHumanSuperbioMetabolism(molecule);
					System.out.println(bts.size() + " biotransformations");
					biotransformations.addAll(bts); 	
				}
				catch(Exception e){
					System.err.println("Could not predicted metabolism for molecule nr. " + nr);
					System.err.println(e.getLocalizedMessage());
				}
			}	
		}
//		this.ecb.saveBioTransformationsToSDF(biotransformations, outputFileName);
		this.ecb.saveBioTransformationProductsToSdf(Utilities.selectUniqueBiotransformations(biotransformations), outputFileName, this.combinedReactionsHash, annotate);
	}	
	
	
	public void simulateHumanSuperbioMetabolismFromSDFtoSingleCSV(String inputFileName, String outputFileName, boolean annotate) throws Exception {
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
				
		int nr = 0;
		IAtomContainerSet containers = FileUtilities.parseSdf(inputFileName);
		if(!containers.isEmpty()){
			for(IAtomContainer molecule : containers.atomContainers()){
				nr++;
				System.out.println("molecule nr. " +nr);
				
				try{
					String identifier = molecule.getProperty(CDKConstants.TITLE);
					if(identifier == null){
						identifier = molecule.getProperty("Name");
						if(identifier == null){
							identifier = molecule.getProperty("$MolName"); 
							if(identifier == null){
								identifier = molecule.getProperty("InChiKey");
								if(identifier == null){
									identifier = this.ecb.inchiGenFactory.getInChIGenerator(molecule).getInchiKey();
								}
							}

						}
					}
					
					ArrayList<Biotransformation> bts = this.simulateHumanSuperbioMetabolism(molecule);
					System.out.println(bts.size() + " biotransformations");
					biotransformations.addAll(bts); 	
				}
				catch(Exception e){
					System.err.println("Could not predicted metabolism for molecule nr. " + nr);
					System.err.println(e.getLocalizedMessage());
				}
			}	
		}
//		this.ecb.saveBioTransformationsToSDF(biotransformations, outputFileName);
		this.ecb.saveBioTransformationProductsToCSV(Utilities.selectUniqueBiotransformations(biotransformations), outputFileName, this.combinedReactionsHash, annotate);
	}
	
	
	public ArrayList<Biotransformation> simulateOneStepAllHuman(IAtomContainer target) throws Exception{
		return simulateOneStepAllHuman(target, 0.5);
	}
	
	public ArrayList<Biotransformation> simulateOneStepAllHuman(IAtomContainer target, double scoreThreshold) throws Exception {
		 ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		 IAtomContainer molecule = ChemStructureManipulator.standardizeMoleculeWithCopy(target);
//		 AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
//		 System.out.println("SMILES BEFORE STANDARDIZATION: " + this.smiGen.create(target));		 
//		 System.out.println("SMILES AFTER STANDARDIZATION: " + this.smiGen.create(molecule));
//
//		 System.out.println("INCHIKEY BEFORE STANDARDIZATION: " + this.ecb.inchiGenFactory.getInChIGenerator(target).getInchiKey());		 
//		 System.out.println("INCHIKEY AFTER STANDARDIZATION: " +  this.ecb.inchiGenFactory.getInChIGenerator(molecule).getInchiKey());
		 
		 IAtomContainerSet products = DefaultChemObjectBuilder
					.getInstance().newInstance(IAtomContainerSet.class);
		 
		
//		System.out.println("Is Biotransformer valid: " + ChemStructureExplorer.isBioTransformerValid(target));
		
		if(ChemStructureExplorer.isBioTransformerValid(molecule)){
			
			ChemicalClassName clname = ChemicalClassFinder.findChemicalClass(molecule);
			
			
			
			if( !(clname == ChemicalClassName.ETHER_LIPID || clname == ChemicalClassName.GLYCEROLIPID || 
					clname == ChemicalClassName.GLYCEROPHOSPHOLIPID ||
					clname == ChemicalClassName.SPHINGOLIPID || clname == ChemicalClassName.GLYCEROL_3_PHOSPHATE_INOSITOL ) ) {
				
				/*
				 *  Apply ECBased
				 *  
				 */	
//				System.out.println("Predict ECBased");
				biotransformations.addAll(this.ecb.simulateECBasedPhaseIMetabolismChain(molecule, true, true, 1, scoreThreshold));
//				biotransformations.addAll(this.ecb.applyEcBasedDeconjuations(molecule, true, true, 3));
//				biotransformations.addAll(this.ecb.applyEcBasedTransformations(molecule, true, true, 0.5));
//				biotransformations.addAll(this.ecb.applyEcBasedConjugations(molecule, true, true, 0.5));
				
				/*
				 *  Apply CYp450
				 *  
				 */
//				System.out.println("Predict CYP450");
				biotransformations.addAll(this.cyb.predictCyp450Biotransformations(molecule, true, true, scoreThreshold));
							
				/*
				 * Human gut metabolism
				 */
//				System.out.println("Predict Human gut metabolism");
				biotransformations.addAll(this.hgb.simulateGutMicrobialMetabolism(molecule, true, true, 1, scoreThreshold));
				
				/*
				 *  Apply Phase II metabolism
				 */	
//				System.out.println("Predict Human Phase II metabolism");
				products.addAtomContainer(molecule);
				biotransformations.addAll(this.p2b.applyPhase2TransformationsChainAndReturnBiotransformations(products, true, true, true, 1, scoreThreshold));
				
				
			
			} else {
				if(ChemicalClassFinder.isEtherLipid(molecule)){
					biotransformations.addAll(this.ecb.applyPathwaySpecificBiotransformations(molecule, MPathwayName.ETHER_LIPID_METABOLISM, true, true, scoreThreshold));
				}
				if(ChemicalClassFinder.isGlyceroLipid(molecule)){
					biotransformations.addAll(this.ecb.applyPathwaySpecificBiotransformations(molecule, MPathwayName.GLYCEROLIPID_METABOLISM, true, true, scoreThreshold));
				}
				if(ChemicalClassFinder.isGlycerol_3_PhosphateInositol(molecule)){
					biotransformations.addAll(this.ecb.applyPathwaySpecificBiotransformations(molecule, MPathwayName.INOSITOL_PHOSPHATE_METABOLISM, true, true, scoreThreshold));
				}								
				if(ChemicalClassFinder.isGlycerophosphoLipid(molecule)){
					biotransformations.addAll(this.ecb.applyPathwaySpecificBiotransformations(molecule, MPathwayName.GLYCEROPHOSPHOLIPID_METABOLISM, true, true, scoreThreshold));
				}		
				if(ChemicalClassFinder.isSphingoLipid(molecule)){
//					System.out.println("Is Sphingolipid");
					biotransformations.addAll(this.ecb.applyPathwaySpecificBiotransformations(molecule, MPathwayName.SPHINGOLIPID_METABOLISM, true, true, scoreThreshold));
				}
				if(ChemicalClassFinder.isC24BileAcid(molecule) || ChemicalClassFinder.isC23BileAcid(molecule)){
//					System.out.println("Is Sphingolipid");
					biotransformations.addAll(this.ecb.applyPathwaySpecificBiotransformations(molecule, MPathwayName.BILE_ACID_METABOLISM, true, true, scoreThreshold));
				}
				
			}			
		}


		return Utilities.selectUniqueBiotransformations(biotransformations);
	}
	
	public ArrayList<Biotransformation> applyPathwaySpecificBiotransformations(IAtomContainer molecule, MPathwayName pathway, boolean preprocess, boolean filter, double scoreThreshold) throws Exception{
		return this.ecb.applyPathwaySpecificBiotransformations(molecule, pathway, preprocess, filter, scoreThreshold);
	}
//	public ArrayList<Biotransformation> simulateECBasedMetabolismChain(IAtomContainer molecule, int nrOfSteps, double scoreThreshold) throws Exception{
//		return this.hgb.simulateGutMicrobialMetabolism(molecule, true, true, nrOfSteps, scoreThreshold);
//	}
	
	public ArrayList<Biotransformation> simulateGutMicrobialMetabolism(IAtomContainer molecule, int nrOfSteps, double scoreThreshold) throws Exception{
		return this.hgb.simulateGutMicrobialMetabolism(molecule, true, true, nrOfSteps, scoreThreshold);
	}
	
	public ArrayList<Biotransformation> applyPhaseIITransformationsChainAndReturnBiotransformations(IAtomContainerSet products, int nrOfSteps, double scoreThreshold) throws Exception{
		return this.p2b.applyPhase2TransformationsChainAndReturnBiotransformations(products, true, true, true, nrOfSteps, scoreThreshold);
	}
	
	public ArrayList<Biotransformation> simulateECBasedPhaseIMetabolismChain(IAtomContainerSet molecules, boolean preprocess, boolean filter,int nrOfSteps, Double scoreThreshold) throws Exception{
		return this.ecb.simulateECBasedPhaseIMetabolismChain(molecules, true, true, nrOfSteps, scoreThreshold);
	}	
	
//	public void simulateHumanMetabolismOneStepFromSDFtoSingleSDF(String inputFileName, String outputFileName) throws Exception{
//		simulateHumanMetabolismOneStepFromSDFtoSingleSDF(inputFileName, outputFileName, 0.5);
//	}	
//	
//	public void simulateHumanMetabolismOneStepFromSDFtoSingleSDF(String inputFileName, String outputFileName, double screThreshold) throws Exception{
//		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
//		
//		int nr = 0;
//		IAtomContainerSet containers = FileUtils.parseSdf(inputFileName);
//		for(IAtomContainer atc : containers.atomContainers()){
//			biotransformations.addAll(this.simulateOneStepAllHuman(atc, screThreshold));
//		}
//	
//		this.ecb.saveBioTransformationsToSDF(biotransformations, outputFileName);
//	}
	
	
	public ArrayList<Biotransformation> simulateOneStepAllHuman(IAtomContainerSet targets, double scoreThreshold) throws Exception {
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		for(IAtomContainer target : targets.atomContainers()){
			biotransformations.addAll(simulateOneStepAllHuman(target, scoreThreshold));
		}
		
		return Utilities.selectUniqueBiotransformations(biotransformations);
	}
	
	public ArrayList<Biotransformation> predictAllHumanBiotransformationChain(IAtomContainer substrate, int nrOfSteps, double threshold) throws Exception{
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		IAtomContainerSet containers = DefaultChemObjectBuilder
				.getInstance().newInstance(IAtomContainerSet.class);
		
		containers.addAtomContainer(substrate);
		int counter = 0;
		while(nrOfSteps>0){
			
			counter++;
			System.out.println("\nStep: " + counter + "\n");
			ArrayList<Biotransformation> currentBiotransformations = simulateOneStepAllHuman(containers, threshold);
			nrOfSteps--;
			if(!currentBiotransformations.isEmpty()){
				biotransformations.addAll(currentBiotransformations);
				containers.removeAllAtomContainers();
				containers = this.ecb.extractAtomContainer(currentBiotransformations);				
			}
			else{
				break;
			}
		}

		return Utilities.selectUniqueBiotransformations(biotransformations);
	}	

	
	public void predictMetabolismAllHumanFromSDFtoSDF(String inputFileName, String outputFileName, boolean annotate) throws Exception{
		predictMetabolismAllHumanFromSDFAndSavetoSDF(inputFileName, outputFileName, 1, 0.5, annotate);
	}

		
	public void predictAllHumanBiotransformationChainAndSaveToSDF(IAtomContainer substrate, int nrOfSteps, double threshold, String outputFileName, boolean annotate) throws Exception{
		ArrayList<Biotransformation> biotransformations= predictAllHumanBiotransformationChain(substrate, nrOfSteps, threshold);
		this.ecb.saveBioTransformationProductsToSdf(Utilities.selectUniqueBiotransformations(biotransformations), outputFileName, this.combinedReactionsHash, annotate);
	}
		
	public void predictMetabolismAllHumanFromSDFAndSavetoSDF(String inputFileName, String outputFileName, int nrOfSteps, boolean annotate) throws Exception{
			predictMetabolismAllHumanFromSDFAndSavetoSDF(inputFileName, outputFileName, nrOfSteps, 0.5, annotate);
//			ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
//			
//			int nr = 0;
//			IAtomContainerSet containers = FileUtils.parseSdf(inputFileName);
//			for(IAtomContainer atc : containers.atomContainers()){
//				biotransformations.addAll(this.predictAllHumanBiotransformationChain(atc, nrOfSteps, 0.5));
//			}	
//			this.ecb.saveBioTransformationsToSDF(biotransformations, outputFileName);
	}	
		
	public void predictMetabolismAllHumanFromSDFAndSavetoSDF(String inputFileName, String outputF, double screThreshold, boolean annotate) throws Exception{
		predictMetabolismAllHumanFromSDFAndSavetoSDF(inputFileName, outputF, 1, screThreshold, annotate);
//		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
//		
//		int nr = 0;
//		IAtomContainerSet containers = FileUtils.parseSdf(inputFileName);
//		for(IAtomContainer atc : containers.atomContainers()){
//			biotransformations.addAll(this.predictAllHumanBiotransformationChain(atc, 1, screThreshold));
//		}	
//		this.ecb.saveBioTransformationsToSDF(biotransformations, outputFileName);
	}	


	public void predictMetabolismAllHumanFromSDFAndSavetoSDF(String inputFileName, String outputFolder, int nrOfSteps, double scoreThreshold, boolean annotate) throws Exception{
			
		int nr = 0;
		IAtomContainerSet containers = FileUtilities.parseSdf(inputFileName);
		for(IAtomContainer atc : containers.atomContainers()){
			try{
				String identifier = atc.getProperty(CDKConstants.TITLE);
				if(identifier == null){
					identifier = atc.getProperty("Name");
					if(identifier == null){
						identifier = atc.getProperty("$MolName"); 
						if(identifier == null){
							identifier = atc.getProperty("InChiKey");
							if(identifier == null){
								identifier = this.ecb.inchiGenFactory.getInChIGenerator(atc).getInchiKey();
							}
						}

					}
				}
				identifier = identifier.replace(":", "-").replace("/", "_");
				ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
				biotransformations.addAll(this.predictAllHumanBiotransformationChain(atc, nrOfSteps, scoreThreshold));
				this.ecb.saveBioTransformationProductsToSdf(Utilities.selectUniqueBiotransformations(biotransformations), outputFolder + "/" + identifier + "_BioT_allHuman_metabolites.sdf", this.combinedReactionsHash, annotate);
			}
			catch(Exception e){
				System.err.println("Could not predicted metabolism for molecule nr. " + nr);
				System.err.println(e.getLocalizedMessage());
			}
		}		
	}	

	
	public void predictMetabolismAllHumanFromSDFAndSavetoSingleSDF(String inputFileName, String outputFileName, int nrOfSteps, double scoreThreshold, boolean annotate) throws Exception{
		
		int nr = 0;
		IAtomContainerSet containers = FileUtilities.parseSdf(inputFileName);
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		for(IAtomContainer atc : containers.atomContainers()){
			try{
				String identifier = atc.getProperty(CDKConstants.TITLE);
				if(identifier == null){
					identifier = atc.getProperty("Name");
					if(identifier == null){
						identifier = atc.getProperty("$MolName"); 
						if(identifier == null){
							identifier = atc.getProperty("InChiKey");
							if(identifier == null){
								identifier = this.ecb.inchiGenFactory.getInChIGenerator(atc).getInchiKey();
							}
						}

					}
				}
				identifier = identifier.replace(":", "-").replace("/", "_");
				
				biotransformations.addAll(this.predictAllHumanBiotransformationChain(atc, nrOfSteps, scoreThreshold));
			}
				
			catch(Exception e){
				System.err.println("Could not predicted metabolism for molecule nr. " + nr);
				System.err.println(e.getLocalizedMessage());
			}
		}
		this.ecb.saveBioTransformationProductsToSdf(Utilities.selectUniqueBiotransformations(biotransformations), outputFileName, this.combinedReactionsHash, annotate);
		
	}	

	public IAtomContainer standardizeMoleculeWithCopy(IAtomContainer target) throws Exception{
		return  ChemStructureManipulator.standardizeMoleculeWithCopy(target);
	}
	

	public ArrayList<Biotransformation> simulateGutMicrobialMetabolismHydrolysisAndReduction(IAtomContainerSet molecules, boolean preprocess, boolean filter , int nrOfSteps, Double scoreThreshold) throws Exception{
		return this.hgb.simulateGutMicrobialMetabolismHydrolysisAndReduction(molecules, preprocess, filter , nrOfSteps, scoreThreshold);
	}
	
	public IAtomContainerSet extractAtomContainer(ArrayList<Biotransformation> biotransformations) throws Exception{
		return this.ecb.extractAtomContainer(biotransformations);
	}
	

	public IAtomContainerSet extractAtomContainerWithTransformationData(ArrayList<Biotransformation> biotransformations, boolean annotate) throws Exception{
		return this.ecb.extractAtomContainerWithTransformationData(biotransformations, this.combinedReactionsHash, annotate);
	}
	
	public void saveBioTransformationProductsToSdf(ArrayList<Biotransformation> biotransformations, String outputFileName, boolean annotate) throws Exception {
		try{
			this.ecb.saveBioTransformationProductsToSdf(Utilities.selectUniqueBiotransformations(biotransformations), outputFileName, this.combinedReactionsHash, annotate);
		}
		catch(Exception e){
			System.err.println(e.getLocalizedMessage());
		}		
	}
	public void saveBioTransformationProductsToCSV(ArrayList<Biotransformation> biotransformations, String outputFileName, boolean annotate) throws Exception {
		try{
			this.ecb.saveBioTransformationProductsToCSV(Utilities.selectUniqueBiotransformations(biotransformations), outputFileName, this.combinedReactionsHash, annotate);
		}
		catch(Exception e){
			System.err.println(e.getLocalizedMessage());
		}		
	}	
	
	
	
//	public IAtomContainerSet extractAtomContainerWithTransformationData(ArrayList<Biotransformation> biotransformations) throws Exception{
////		AtomContainerSet acontainers = new AtomContainerSet();
//		IAtomContainerSet acontainers = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
//		LinkedHashMap<String, IAtomContainer> hMap = new LinkedHashMap<String, IAtomContainer>();
//		int metaboliteID = 0;
//		for(Biotransformation b : biotransformations){
//			for(IAtomContainer ac : b.getProducts().atomContainers()){
//				IAtomContainer hash_ac;				
//				String ikey = ac.getProperty("InChIKey");
//				if(hMap.containsKey(ikey)) {
//					hash_ac = hMap.get(ikey);
//					Set<String> r = new HashSet<String>(Arrays.asList(hash_ac.getProperty("Reactions").toString().split("\n")));
////					System.out.println("reactions: " +  r);
////					if(b.getEnzymeNames() == null){
////						System.out.println("b.getEnzymeNames() == null");
//						r.add(b.getReactionType().toString());
////					}
////					else 
//						if(b.getEnzymeNames().size()>0){
////						System.out.println("b.getEnzymeNames().size() > 0 ");
////						r.add(b.getReactionType().toString() + " (" + StringUtils.join(b.getEnzymeNames(), ", ") + ")");
//						ac.setProperty("Enzyme(s)", StringUtils.join(b.getEnzymeNames(), "\n"));
//						
////						System.out.println("Molecule: " + ac.getProperty("InChI"));
////						System.out.println("Reaction: " + ac.getProperty("Reactions"));
////						System.out.println("Enzymes: " + ac.getProperty("Enzyme(s)"));
//						
//					}
//
////					hash_ac.setProperty("Reactions", StringUtils.join(r, "\n"));	
////					Set<String> e = new HashSet<String>(Arrays.asList(hash_ac.getProperty("Enzymes").toString().split("\n")));
////					ArrayList al =new ArrayList<String>(e);
////					hash_ac.setProperty("Enzymes", StringUtils.join(al, "\n"));
//					
////					e.add(b.getReactionType().toString());
//					
//				}
//				else {
//					ac.setProperty("Molecular formula", ChemStructureExplorer.getMolecularFormula(ac));
////					ac.setProperty("Reactions", b.getReactionType().toString());
////					System.out.println(b.getReactionType().toString());
////					System.err.println(this.reactionsHash);
////					System.err.println(this.reactionsHash.size());
////					System.out.println(this.reactionsHash.get(b.getReactionType().toString()).getComonName());
////					System.out.println(this.reactionsHash.get(b.getReactionType().toString()).getBTRMID());
//					ac.setProperty("Reactions", this.reactionsHash.get(b.getReactionType()).getComonName() + 
//							" (" + this.reactionsHash.get(b.getReactionType()).getBTRMID() +")" );
//					if(b.getEnzymeNames().size()>0){
//						ac.setProperty("Enzyme(s)", StringUtils.join(b.getEnzymeNames(),"\n"));
//						
////						System.out.println("Molecule: " + ac.getProperty("InChI"));
////						System.out.println("Reaction 2: " + ac.getProperty("Reactions"));
////						System.out.println("Enzymes 2: " + ac.getProperty("Enzyme(s)"));
//					}
//					metaboliteID++;
////					ac.setProperty("Metabolite ID", "BTM" + String.format("%04d", metaboliteID));
//					ac.setProperty(CDKConstants.TITLE, "BTM" + String.format("%05d", metaboliteID));
//					hMap.put(ikey, ac);
//					
//				}
//				
//				
//			if(b.getSubstrates().getAtomContainerCount() == 1){
//
//				
//				String tt = (String) b.getSubstrates().getAtomContainer(0).getProperty(CDKConstants.TITLE);
//				if(tt == null){
//					tt = (String) b.getSubstrates().getAtomContainer(0).getProperty("Name");
//					if(tt == null){
//						tt = (String) b.getSubstrates().getAtomContainer(0).getProperty("Metabolite ID");
//					}
//				}
//				
//				ChemStructureExplorer.addPhysicoChemicalProperties(ac);
//				
//				if(b.getSubstrates().getAtomContainer(0).getProperty("Major Isotope Mass") == null){
//					ChemStructureExplorer.addPhysicoChemicalProperties(b.getSubstrates().getAtomContainer(0));
//				}
//				
//				
//				ac.setProperty("BioSystem", b.getBioSystemName());
////				ac.setProperty("Precursor", b.getSubstrates().getAtomContainer(0).getProperty(tt));
//				
//				
//				if(hMap.containsKey(b.getSubstrates().getAtomContainer(0).getProperty("InChIKey"))){
//					ac.setProperty("Precursor ID", tt);
//				}
////				else{
////					
////				}
////				if(b.getSubstrates().getAtomContainer(0).getProperty("Major Isotope Mass")){
////					ac.setProperty("Precursor ID", b.getSubstrates().getAtomContainer(0).getProperty("Metabolite ID"));
////				}
//				
//				
//				ac.setProperty("Precursor InChI", b.getSubstrates().getAtomContainer(0).getProperty("InChI"));
//				ac.setProperty("Precursor InChIKey", b.getSubstrates().getAtomContainer(0).getProperty("InChIKey"));
////				ac.setProperty("Precursor XLogP", b.getSubstrates().getAtomContainer(0).getProperty("XLogP"));
//				ac.setProperty("Precursor ALogP", b.getSubstrates().getAtomContainer(0).getProperty("ALogP"));
//				ac.setProperty("Precursor Major Isotope Mass", b.getSubstrates().getAtomContainer(0).getProperty("Major Isotope Mass"));
////				ac.setProperty("Precursor Molecular weight", b.getSubstrates().getAtomContainer(0).getProperty("Molecular weight"));
//				
//			}
//			}
//		}
//		ArrayList<IAtomContainer> at =  new ArrayList<IAtomContainer>( hMap.values());
////		System.err.println(at.get(0).getClass());
//		for(IAtomContainer a : at){
//			acontainers.addAtomContainer(a);
//		}
//		
//		return ChemStructureExplorer.uniquefy(acontainers);
//	}
//	

	

	}
