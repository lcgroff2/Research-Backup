/**
 * This class implements the class of human gut biotransformers, which simulate the transformation of molecules
 * by enzymes from microbial species found in the human gut. It implements rules and constraints extracted from or designed upon
 * (1) mining the scientific literature, (2) expert collaboration, and/or (3) experimental validation.
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */


package biotransformer.btransformers;

import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.json.simple.parser.ParseException;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import ambit2.smarts.SMIRKSManager;
import ambit2.smarts.query.SMARTSException;
import ambit2.smarts.query.SmartsPatternCDK;
import biotransformer.biomolecule.Enzyme;
import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.btransformers.Biotransformer;
import biotransformer.transformation.Biotransformation;
import biotransformer.transformation.MReactionSets;
import biotransformer.transformation.MetabolicReaction;
import biotransformer.utils.ChemStructureExplorer;
import biotransformer.utils.ChemStructureManipulator;
import biotransformer.utils.Utilities;

public class HGutBTransformer extends Biotransformer {

	public HGutBTransformer() throws IOException, CDKException{
		super(BioSystemName.GUTMICRO);
		setReactionsList();
	}
	
	/**
	 * Collect the list of metabolic reactions associated with the current biosystem, inferred from the list of enzymes.
	 * The deconjugation enzymes include esterases, alpha-rhamnosidases, beta-glucuronidases, beta-glycosidases.
	 */
	private void setReactionsList(){
		//		this.reactionsList.put("standardizationReactions", MReactionSets.standardizationReactions);
		ArrayList<MetabolicReaction> reductionReact = new ArrayList<MetabolicReaction>() ;
		ArrayList<MetabolicReaction> phaseIIReact = new ArrayList<MetabolicReaction>();
		ArrayList<MetabolicReaction> deconjugationReact = new ArrayList<MetabolicReaction>();
//		this.enzymesByreactionGroups.put("deconjugationReactions", new ArrayList<Enzyme>());
//		this.enzymesByreactionGroups.put("gutMicroReductiveReactions", new ArrayList<Enzyme>());
//		this.enzymesByreactionGroups.put("gutMicroPhaseIIReactions", new ArrayList<Enzyme>());

		ArrayList<Enzyme> reductionEnz = new ArrayList<Enzyme>() ;
		ArrayList<Enzyme> phaseIIEnz = new ArrayList<Enzyme>();
		ArrayList<Enzyme> deconjugationEnz = new ArrayList<Enzyme>();
		
				
		for( Enzyme enz : this.bSystem.getEnzymesList()){
			if(enz.getName().contentEquals("EC_2_1_1_6") || enz.getName().contentEquals("SULFOTRANSFERASE") ||
					enz.getName().contentEquals("UDP_GLUCURONOSYLTRANSFERASE")){
				phaseIIReact.addAll(enz.getReactionSet());
				phaseIIEnz.add(enz);	
//				this.enzymesByreactionGroups.get("gutMicroPhaseIIReactions").add(enz);
			}
			
			else if (
					
					// carboxylesterase
					enz.getName().contentEquals("EC_3_1_1_1") || enz.getName().contentEquals("EC_3_1_1_2") || enz.getName().contentEquals("EC_3_1_1_20") || 
					
					// sulfatase
					
					enz.getName().contentEquals("EC_3_1_6_1") || 
					
					// glycosidase
					enz.getName().contentEquals("EC_3_2_1_X") || enz.getName().contentEquals("EC_3_2_1_20") || 
					enz.getName().contentEquals("EC_3_2_1_21") || enz.getName().contentEquals("EC_3_2_1_23") || 
					enz.getName().contentEquals("EC_3_2_1_31") || enz.getName().contentEquals("EC_3_2_1_40") ||
					enz.getName().contentEquals("EC_3_2_1_147") || enz.getName().contentEquals("FLAVONOID_C_GLYCOSIDASE") || 
					enz.getName().contentEquals("EC_3_1_1_73") ||
					
					// hydroxylase
					enz.getName().contentEquals("EC_3_5_1_24")

									
					
					){
				deconjugationReact.addAll(enz.getReactionSet());
				deconjugationEnz.add(enz);
//				this.enzymesByreactionGroups.get("deconjugationReactions").add(enz);
			}
			
			else{
				reductionReact.addAll(enz.getReactionSet());
//				System.err.println(enz.getName());
				reductionEnz.add(enz);
				
			}
		}
		
		
		
		this.reactionsByGroups.put("gutMicroReductiveReactions", reductionReact);
		this.reactionsByGroups.put("gutMicroPhaseIIReactions", phaseIIReact);
		// Phase O consists of deglycosyltation, de-esterification, desulfation 
		this.reactionsByGroups.put("deconjugationReactions", deconjugationReact );
		
		
//		System.err.println("deconjugationReactions: " + deconjugationEnz.size());
//		System.err.println("gutMicroReductiveReactions: " + reductionEnz.size());
//		System.err.println("gutMicroPhaseIIReactions: " + phaseIIEnz.size());
		this.enzymesByreactionGroups.put("deconjugationReactions", deconjugationEnz);
		this.enzymesByreactionGroups.put("gutMicroReductiveReactions", reductionEnz);
		this.enzymesByreactionGroups.put("gutMicroPhaseIIReactions", phaseIIEnz);
		
//		put("gutMicroReductiveEnzymes", reductionEnz);
//		put("gutMicroPhaseIIEnzymes", phaseIIEnz);
//		// Phase O consists of deglycosyltation, de-esterification, desulfation 
//		put("deconjugationEnzymes", deconjugationEnz );
		
		for(MetabolicReaction m : this.reactionsByGroups.get("gutMicroReductiveReactions")){
			this.reactionsHash.put(m.name, m);
		}
		for(MetabolicReaction m : this.reactionsByGroups.get("gutMicroPhaseIIReactions")){
			this.reactionsHash.put(m.name, m);
		}
		for(MetabolicReaction m : this.reactionsByGroups.get("deconjugationReactions")){
			this.reactionsHash.put(m.name, m);
		}

	}
	/**
	 * returns a linked hash map with the reactions associated with the human gut, in addition to standardization reactions.
	 */
	public LinkedHashMap<String, ArrayList<MetabolicReaction>> getReactionsList(){
		return this.reactionsByGroups;
	}
	
	/**
	 * 
	 * @param target - The molecule to transform
	 * @param preprocess - specify whether to perform molecule preprocessing
	 * @param filter - apply reaction filtering
	 * @return an arraylist of biotransformations, which are instances of the human gut microbial reactions applied to the target, with a threshold of 0.0
	 * @throws Exception - throw any exception
	 */
	public ArrayList<Biotransformation> applyGutMicrobialDeconjugations(IAtomContainer target, boolean preprocess, boolean filter)
			throws Exception {
		return this.applyGutMicrobialDeconjugations(target, preprocess, filter, 0.0);
	}
	
	/**
	 * 
	 * @param target - The molecule to transform
	 * @param preprocess - specify whether to perform molecule preprocessing
	 * @param filter - apply reaction filtering
	 * @return an arraylist of biotransformations, which are instances of the human gut microbial reactions applied to the target, with the set minimum threshold
	 * @throws Exception -  throw any exception
	 */
	public ArrayList<Biotransformation> applyGutMicrobialDeconjugations(IAtomContainer target, boolean preprocess, boolean filter, double scoreThreshold)
			throws Exception {
		
		

		
		if(ChemStructureExplorer.isBioTransformerValid(target)) {
			
//			return this.applyReactionsAndReturnBiotransformations(target, this.reactionsByGroups.get("deconjugationReactions"), preprocess, filter, scoreThreshold);
			return this.metabolizeWithEnzymes(target, this.enzymesByreactionGroups.get("deconjugationReactions"), preprocess, filter, scoreThreshold);
		
		}
		else if(ChemStructureExplorer.isMixture(target) || ChemStructureExplorer.isCompoundInorganic(target)){
//			throw new IllegalArgumentException("The substrate must not be a mixture.");
			throw new IllegalArgumentException("\n\n INVALID COMPOUND:\nThe compound is not valid for BioTransformer. Make sure that the compound is: 1) organic; and 2) not a mixture.");
		}
//		
		else{
			return null;
		}
	}

	/**
	 * 
	 * @param target - The molecule to transform
	 * @param preprocess - specify whether to perform molecule preprocessing
	 * @param filter - apply reaction filtering
	 * @param nr_of_steps -  number of steps
	 * @return an arraylist of biotransformations obtained after the specified number of steps (nr_of_steps), which are 
	 * instances of the human gut microbial metabolic reactions applied to the target, with the minimum threshold of 0.0
	 * @throws Exception - throw any exception
	 */
	public ArrayList<Biotransformation> applyGutMicrobialDeconjugationsChain(IAtomContainer target,
			boolean preprocess, boolean filter, int nr_of_steps) throws Exception{
		return  this.applyGutMicrobialDeconjugationsChain(target, preprocess, filter, nr_of_steps, 0.0);
	}
	
	/**
	 * 
	 * @param target - The molecule to transform
	 * @param preprocess - specify whether to perform molecule preprocessing
	 * @param filter - apply reaction filtering
	 * @param nr_of_steps -  number of steps
	 * @return an arraylist of biotransformations obtained after the specified number of steps (nr_of_steps), which are 
	 * instances of the human gut microbial metabolic reactions applied to the target, with the set minimum threshold
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> applyGutMicrobialDeconjugationsChain(IAtomContainer target,
			boolean preprocess, boolean filter, int nr_of_steps, Double scoreThreshold) throws Exception{
		
//		if(ChemStructureExplorer.isMixture(target) || ChemStructureExplorer.isCompoundInorganic(target)){
////			throw new IllegalArgumentException("The substrate must not be a mixture.");
//			System.err.println("\n\n INVALID COMPOUND:\nThe compound is not valid for BioTransformer. Make sure that the compound is: 1) organic; and 2) not a mixture.");
////			return biotransformations;
//		}
//		
//		if(ChemStructureExplorer.isBioTransformerValid(target)) {
//			return this.applyReactionsChainAndReturnBiotransformations(target, this.reactionsByGroups.get("deconjugationReactions"), preprocess, filter, nr_of_steps, scoreThreshold);
			return this.metabolizeWithEnzymesBreadthFirst(target, this.enzymesByreactionGroups.get("deconjugationReactions"), preprocess, filter, nr_of_steps, scoreThreshold);
//		}else{
//			throw new IllegalArgumentException("The target must not be a mixture.");
//		}
	}
	
	
	
	/**
	 * 
	 * @param target - The molecule to transform
	 * @param preprocess - specify whether to perform molecule preprocessing
	 * @param filter - apply reaction filtering
	 * @return an arraylist of biotransformations, which are instances of the human gut microbial reactions applied to the target, with a threshold of 0.0
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> applyGutMicrobialReductions(IAtomContainer target, boolean preprocess, boolean filter)
			throws Exception {
		return this.applyGutMicrobialReductions(target, preprocess, filter, 0.0);
	}
	
	/**
	 * 
	 * @param target - The molecule to transform
	 * @param preprocess - specify whether to perform molecule preprocessing
	 * @param filter - apply reaction filtering
	 * @param scoreThreshold - minimum threshold for reaction scores
	 * @return an arraylist of biotransformations, which are instances of the human gut microbial reactions applied to the target, with the set minimum threshold
	 * @throws Exception - throw any exception
	 */
	public ArrayList<Biotransformation> applyGutMicrobialReductions(IAtomContainer target, boolean preprocess, boolean filter, double scoreThreshold)
			throws Exception {
		
		try {
			if(ChemStructureExplorer.isCompoundInorganic(target) || ChemStructureExplorer.isMixture(target)){
				throw new IllegalArgumentException(target.getProperty("InChIKey")+ "\nThe substrate must be: 1) organic, and; 2) not a mixture.");
			} else if(ChemStructureExplorer.isBioTransformerValid(target)){
				return this.metabolizeWithEnzymes(target, this.enzymesByreactionGroups.get("gutMicroReductiveReactions"), preprocess, filter, scoreThreshold);
	//			return this.applyReactionsAndReturnBiotransformations(target, this.reactionsByGroups.get("gutMicroReductiveReactions"), preprocess, filter, scoreThreshold);
			}else{
				return null;
			}
		}
		catch(Exception e){
			e.printStackTrace();
			return null;
		}
	}

	/**
	 * 
	 * @param target - The molecule to transform
	 * @param preprocess - specify whether to perform molecule preprocessing
	 * @param filter - apply reaction filtering
	 * @param nr_of_steps -  number of steps
	 * @return an arraylist of biotransformations obtained after the specified number of steps (nr_of_steps), which are 
	 * instances of the human gut microbial metabolic reactions applied to the target, with the minimum threshold of 0.0
	 * @throws Exception - throw any exception
	 */
	public ArrayList<Biotransformation> applyGutMicrobialReductionsChain(IAtomContainer target,
			boolean preprocess, boolean filter, int nr_of_steps) throws Exception{
		return  this.applyGutMicrobialReductionsChain(target, preprocess, filter, nr_of_steps, 0.0);
	}
	
	/**
	 * 
	 * @param target - The molecule to transform
	 * @param preprocess - specify whether to perform molecule preprocessing
	 * @param filter - apply reaction filtering
	 * @param nr_of_steps -  number of steps
	 * @return an arraylist of biotransformations obtained after the specified number of steps (nr_of_steps), which are 
	 * instances of the human gut microbial metabolic reactions applied to the target, with the set minimum threshold
	 * @throws Exception -throw any exception
	 */
	public ArrayList<Biotransformation> applyGutMicrobialReductionsChain(IAtomContainer target,
			boolean preprocess, boolean filter, int nr_of_steps, Double scoreThreshold) throws Exception{
				
		try{
			if(ChemStructureExplorer.isCompoundInorganic(target) || ChemStructureExplorer.isMixture(target)){
				throw new IllegalArgumentException(target.getProperty("InChIKey")+ "\nThe substrate must be: 1) organic, and; 2) not a mixture.");
			} else if(ChemStructureExplorer.isBioTransformerValid(target)) {	
				if( !this.isDeconjugationCandidate(target)){
					return this.metabolizeWithEnzymesBreadthFirst(target, this.enzymesByreactionGroups.get("gutMicroReductiveReactions"), preprocess, filter, nr_of_steps, scoreThreshold);
	//				return this.applyReactionsChainAndReturnBiotransformations(target, this.reactionsByGroups.get("gutMicroReductiveReactions"), preprocess, filter, nr_of_steps, scoreThreshold);			
				}else{
					return null;
				}
			}else{
				return null;
			}
		}
		catch(Exception e){
			e.printStackTrace();
			return null;
		}
		
	}
	
	

	/**
	 * 
	 * @param target - The molecule to transform
	 * @param preprocess - specify whether to perform molecule preprocessing
	 * @param filter - apply reaction filtering
	 * @return an arraylist of biotransformations, which are instances of the human gut microbial reactions applied to the target, with a threshold of 0.0
	 * @throws Exception - throw any exception
	 */
	public ArrayList<Biotransformation> applyGutMicrobialConjugations(IAtomContainer target, boolean preprocess, boolean filter)
			throws Exception {
		return this.applyGutMicrobialConjugations(target, preprocess, filter, 0.0);
	}
	
	/**
	 * 
	 * @param target - The molecule to transform
	 * @param preprocess - specify whether to perform molecule preprocessing
	 * @param filter - apply reaction filtering
	 * @return an arraylist of biotransformations, which are instances of the human gut microbial reactions applied to the target, with the set minimum threshold
	 * @throws Exception - throw any exception
	 */
	public ArrayList<Biotransformation> applyGutMicrobialConjugations(IAtomContainer target, boolean preprocess, boolean filter, double scoreThreshold)
			throws Exception {
		
		try{
			if(ChemStructureExplorer.isCompoundInorganic(target) || ChemStructureExplorer.isMixture(target)){
				throw new IllegalArgumentException(target.getProperty("InChIKey")+ "\nThe substrate must be: 1) organic, and; 2) not a mixture.");
			} 
			else if(ChemStructureExplorer.isBioTransformerValid(target)) {		
				
				ArrayList<Biotransformation> biotransformations  = this.metabolizeWithEnzymes(target, this.enzymesByreactionGroups.get("gutMicroPhaseIIReactions"), preprocess, filter, scoreThreshold);
				
	//			ArrayList<Biotransformation> biotransformations  = this.applyReactionsAndReturnBiotransformations(target, this.reactionsByGroups.get("gutMicroPhaseIIReactions"), preprocess, filter, scoreThreshold);
				
				ArrayList<Biotransformation> selectedBiotransformations = new ArrayList<Biotransformation>();
				
				for(Biotransformation bt : biotransformations) {
					boolean goodBiontransfo =  true;
		
	//				System.out.println("Biotransformation type. " +  bt.getReactionType());
					for( IAtomContainer at : bt.getProducts().atomContainers() ){
						
	//					System.out.println(at.getProperty("InChI"));
	//					at = ChemStructureManipulator.preprocessContainer(at);
						if(ChemStructureExplorer.isInvalidPhase2Metabolite(at)){
	//						System.err.println("is invalid phase 2 metabolite");
							goodBiontransfo = false;					
							break;
						}
					}
					
	//				System.out.println("goodBiontransfo: " + goodBiontransfo);
					if(goodBiontransfo){
						selectedBiotransformations.add(bt);
					}
				}
				
				return selectedBiotransformations;
			} else {
				return null;
			}
		}
			catch(Exception e){
				e.printStackTrace();
				return null;
			}
		
	}

	/**
	 * 
	 * @param target
	 * @param preprocess
	 * @param filter
	 * @param nr_of_steps
	 * @return an arraylist of biotransformations obtained after the specified number of steps (nr_of_steps), which are 
	 * instances of the human gut microbial metabolic reactions applied to the target, with the minimum threshold of 0.0
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> applyGutMicrobialConjugationsChain(IAtomContainer target,
			boolean preprocess, boolean filter, int nr_of_steps) throws Exception{
		try{
			if(ChemStructureExplorer.isCompoundInorganic(target) || ChemStructureExplorer.isMixture(target)){
				throw new IllegalArgumentException(target.getProperty("InChIKey")+ "\nThe substrate must be: 1) organic, and; 2) not a mixture.");
			} 
			else if(ChemStructureExplorer.isBioTransformerValid(target))  {				
				return  this.applyGutMicrobialConjugationsChain(target, preprocess, filter, nr_of_steps, 0.0);
			}
			else {
				return null;
			}
		}
		catch(Exception e){
			e.printStackTrace();
			return null;
		}
	}
	
	/**
	 * 
	 * @param target
	 * @param preprocess
	 * @param filter
	 * @param nr_of_steps
	 * @return an arraylist of biotransformations obtained after the specified number of steps (nr_of_steps), which are 
	 * instances of the human gut microbial metabolic reactions applied to the target, with the set minimum threshold
	 * @throws Exception
	 */
	public ArrayList<Biotransformation> applyGutMicrobialConjugationsChain(IAtomContainer target,
			boolean preprocess, boolean filter, int nr_of_steps, Double scoreThreshold) throws Exception{
	
		ArrayList<Biotransformation> biotransformations  = this.metabolizeWithEnzymesBreadthFirst(target, this.enzymesByreactionGroups.get("gutMicroPhaseIIReactions"), preprocess, filter, nr_of_steps, scoreThreshold);
//		ArrayList<Biotransformation> biotransformations  = this.applyReactionsChainAndReturnBiotransformations(target, 
//				this.reactionsByGroups.get("gutMicroPhaseIIReactions"), preprocess, filter, nr_of_steps, scoreThreshold);
		ArrayList<Biotransformation> selectedBiotransformations = new ArrayList<Biotransformation>();
		
		for(Biotransformation bt : biotransformations) {
			boolean goodBiontransfo =  true;

//			System.out.println("Biotransformation type. " +  bt.getReactionType());
			for( IAtomContainer at : bt.getProducts().atomContainers() ){
				
//				System.out.println(at.getProperty("InChI"));
//				AtomContainerManipulator.convertImplicitToExplicitHydrogens(at);
//				at = ChemStructureManipulator.preprocessContainer(at);
				if(ChemStructureExplorer.isInvalidPhase2Metabolite(at)){
//					System.err.println("is invalid phase 2 metabolite");
					goodBiontransfo = false;					
					break;
				}
			}
			
//			System.out.println("goodBiontransfo: " + goodBiontransfo);
			if(goodBiontransfo){
				selectedBiotransformations.add(bt);
			}
		}
		
		
		return selectedBiotransformations;
	
	
	}
	
	
	public ArrayList<Biotransformation> simulateGutMicrobialMetabolismHydrolysisAndReduction(IAtomContainer target,
			boolean preprocess, boolean filter, int nr_of_steps, Double scoreThreshold) throws Exception{
			
		try {
//			ChemStructureExplorer.addInChIandKey(target);
			if(ChemStructureExplorer.isCompoundInorganic(target) || ChemStructureExplorer.isMixture(target)){
				throw new IllegalArgumentException(target.getProperty("InChIKey")+ "\nThe substrate must be: 1) organic, and; 2) not a mixture.");
			} 
			else if(ChemStructureExplorer.isBioTransformerValid(target)){
//				System.out.println("\n\n===========================================");
//				System.out.println("Predicting human gut microbial metabolism");
//				System.out.println("===========================================\n\n");
				ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
//				if(ChemStructureExplorer.isPolyphenolOrDerivative(target)) {
					
					ArrayList<Biotransformation> reductionBiotransformations = new ArrayList<Biotransformation>();
					IAtomContainerSet products = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
					IAtomContainer st = ChemStructureManipulator.standardizeMoleculeWithCopy(target, true);
//					System.out.println(this.smiGen.create(st));
					
			//		biotransformations = this.applyReactionAtOnceAndReturnBiotransformations(target,
			//				this.reactionsList.get("deconjugationReactions"), preprocess, filter, nr_of_steps, scoreThreshold);
						
					
					if(this.isDeconjugationCandidate(st)){
//						System.out.println("Must be deconjugated");
//						biotransformations = this.applyGutMicrobialDeconjugationsChain(target,
//								preprocess, filter, 5, scoreThreshold);
						biotransformations = this.applyGutMicrobialDeconjugationsChain(target,
								preprocess, filter, 1, scoreThreshold);
						products = this.extractAtomContainer(biotransformations);
//						System.out.println("Metabolites after deconjugation: " + products.getAtomContainerCount());
					}
			
					if(products.getAtomContainerCount() == 0){
						// In case the original compound was ready for reduction,e.g. if it did not have any sulfate, glycosyl, etc..
						products.addAtomContainer(target);
					}
			
					for(IAtomContainer s : products.atomContainers()){
						IAtomContainer sc = ChemStructureManipulator.standardizeMoleculeWithCopy(s, true);
		//				 System.out.println("The compound " + this.smiGen.isomeric().create(sc));
		//				 System.out.println("Is a deconjugation candidate " + this.isDeconjugationCandidate(sc));
						
//						if((!this.isDeconjugationCandidate(sc)) && ChemStructureExplorer.isPhaseIPolyphenolCandidateOrDerivative(sc)){
						if((!this.isDeconjugationCandidate(sc))){
//											System.out.println("Is a metabolizable polyphenol\n");
							ArrayList<Biotransformation> acs = applyGutMicrobialReductionsChain(sc,
									preprocess, filter, nr_of_steps, scoreThreshold);
							reductionBiotransformations.addAll(acs);
						}
					}
					
					biotransformations.addAll(reductionBiotransformations);
					IAtomContainerSet reductionProducts = this.extractAtomContainer(reductionBiotransformations);

//				}
	
				
				
				
	//			for(IAtomContainer p : reductionProducts.atomContainers()){
	//				IAtomContainer sp = standardizeMoleculeWithCopy(p, true);
	//				ArrayList<Biotransformation> phaseIIBiotransformations = applyGutMicrobialConjugationsChain(sp,
	//						preprocess, filter, 1, scoreThreshold);
	//				biotransformations.addAll(phaseIIBiotransformations);
	//			}
		
				return Utilities.selectUniqueBiotransformations(reductionBiotransformations);
	
			} else{
				return null;
			}
		} catch (Exception iae) {
			System.err.println(iae.getLocalizedMessage());
//			return biotransformations;
			return null;
		}
	}
	
	
	
	public ArrayList<Biotransformation> simulateGutMicrobialMetabolismHydrolysisAndReduction(IAtomContainerSet targets, boolean preprocess, 
			boolean filter, int nr_of_steps, Double scoreThreshold) throws Exception{
		
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
//		System.out.println("TARGETS");
//		System.out.println(targets == null);
		for(IAtomContainer atc : targets.atomContainers()){
			ArrayList<Biotransformation> bts = this.simulateGutMicrobialMetabolismHydrolysisAndReduction(atc, preprocess, 
					filter, nr_of_steps, scoreThreshold);
			if(bts != null && bts.size()>0){
				biotransformations.addAll(bts);
			}
		}
		
		return biotransformations;
		
	}
	

	public ArrayList<Biotransformation> simulateGutMicrobialMetabolism(IAtomContainer target,
			boolean preprocess, boolean filter, int nr_of_steps, Double scoreThreshold) throws Exception{
		
		try{
			ChemStructureExplorer.addInChIandKey(target);
			if(ChemStructureExplorer.isCompoundInorganic(target) || ChemStructureExplorer.isMixture(target)){
				throw new IllegalArgumentException(target.getProperty("InChIKey")+ "\nThe substrate must be: 1) organic, and; 2) not a mixture.");
			} else if(ChemStructureExplorer.isBioTransformerValid(target)) {		

				ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
				
//				if(ChemStructureExplorer.isPolyphenolOrDerivative(target)){
					ArrayList<Biotransformation> reductionBiotransformations = new ArrayList<Biotransformation>();
					IAtomContainerSet products = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
					IAtomContainer st = ChemStructureManipulator.standardizeMoleculeWithCopy(target, true);
//					System.out.println(this.smiGen.create(st));
					
			//		biotransformations = this.applyReactionAtOnceAndReturnBiotransformations(target,
			//				this.reactionsList.get("deconjugationReactions"), preprocess, filter, nr_of_steps, scoreThreshold);
					
//					System.out.println("Is polyphenol? " + ChemStructureExplorer.isPolyphenolOrDerivative(st));
//					if(ChemStructureExplorer.isPolyphenolOrDerivative(st)){
						if(this.isDeconjugationCandidate(st)){
//							System.out.println("Must be deconjugated");
//							biotransformations = this.applyGutMicrobialDeconjugationsChain(target,
//									preprocess, filter, 5, scoreThreshold);
							biotransformations = this.applyGutMicrobialDeconjugationsChain(target,
									preprocess, filter, 1, scoreThreshold);
							products = this.extractAtomContainer(biotransformations);
//							System.err.println("Metabolites after deconjugation: " + products.getAtomContainerCount());
						}
				
						if(products.getAtomContainerCount() == 0){
							// In case the original compound was ready for reduction,e.g. if it did not have any sulfate, glycosyl, etc..
							products.addAtomContainer(target);
						}
				
						for(IAtomContainer s : products.atomContainers()){
							IAtomContainer sc = ChemStructureManipulator.standardizeMoleculeWithCopy(s, true);
//							 System.out.println("The compound " + this.smiGen.isomeric().create(sc));
//							 System.out.println("Is a deconjugation candidate " + this.isDeconjugationCandidate(sc));
							
//							if((!this.isDeconjugationCandidate(sc)) && ChemStructureExplorer.isPhaseIPolyphenolCandidateOrDerivative(sc)){
							if((!this.isDeconjugationCandidate(sc))){	
							
//												System.out.println("Is a metabolizable polyphenol\n");
								ArrayList<Biotransformation> acs = applyGutMicrobialReductionsChain(sc,
										preprocess, filter, nr_of_steps, scoreThreshold);
								reductionBiotransformations.addAll(acs);
							}
						}
						
						biotransformations.addAll(reductionBiotransformations);
						IAtomContainerSet reductionProducts = this.extractAtomContainer(reductionBiotransformations);
//						System.out.println("Nr. of reduction biotransformations: " + reductionProducts.getAtomContainerCount());
						
						for(IAtomContainer p : reductionProducts.atomContainers()){
							IAtomContainer sp = ChemStructureManipulator.standardizeMoleculeWithCopy(p, true);
							ArrayList<Biotransformation> phaseIIBiotransformations = applyGutMicrobialConjugationsChain(sp,
									preprocess, filter, 1, scoreThreshold);
							biotransformations.addAll(phaseIIBiotransformations);
						}
				
//						System.out.println("Nr. of biotransformations: " + biotransformations.size());				
//					}				
//				}
			

				return Utilities.selectUniqueBiotransformations(biotransformations);

			}
			else{
				return null;
			}
		}catch(IllegalArgumentException e){
			e.printStackTrace();
			return null;
		}
	}
	

	
	public ArrayList<Biotransformation> simulateGutMicrobialMetabolism(IAtomContainerSet targets, boolean preprocess, 
			boolean filter, int nr_of_steps, Double scoreThreshold) throws Exception{
		
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();
		for(IAtomContainer atc : targets.atomContainers()){
			ArrayList<Biotransformation> bts = this.simulateGutMicrobialMetabolism(atc, preprocess, 
					filter, nr_of_steps, scoreThreshold);
			if(bts != null && bts.size()>0){
				biotransformations.addAll(bts);
			}
		}
		
		return Utilities.selectUniqueBiotransformations(biotransformations);
		
	}
	
	public ArrayList<Biotransformation> simulateGutMicrobialMetabolism(String targetsFileNameInSDF, boolean preprocess, 
			boolean filter, int nr_of_steps, Double scoreThreshold) throws Exception{
		
		ArrayList<Biotransformation> biotransformations = new ArrayList<Biotransformation>();			
		IteratingSDFReader sdfr = new IteratingSDFReader(new FileReader(targetsFileNameInSDF), this.builder);
		
		while (sdfr.hasNext()){
			IAtomContainer mol = sdfr.next();
			ArrayList<Biotransformation> bts = this.simulateGutMicrobialMetabolism(mol, preprocess, 
					filter, nr_of_steps, scoreThreshold);
			if(bts!= null && bts.size()>0){
				biotransformations.addAll(bts);
			}
		}
		
		return Utilities.selectUniqueBiotransformations(biotransformations);
	}

	public void simulateGutMicrobialMetabolismAndSave(String targetsFileNameInSDF, boolean preprocess, 
			boolean filter, int nr_of_steps, Double scoreThreshold, String metabolitesFileNameInSDF) throws Exception{

		simulateGutMicrobialMetabolismAndSave(targetsFileNameInSDF, preprocess, 
				filter, nr_of_steps, scoreThreshold, metabolitesFileNameInSDF, false);
	}
	
	
	public void simulateGutMicrobialMetabolismAndSave(String targetsFileNameInSDF, boolean preprocess, 
			boolean filter, int nr_of_steps, Double scoreThreshold, String metabolitesFileNameInSDF, boolean annotate) throws Exception{
		
		ArrayList<Biotransformation> biotransformations = simulateGutMicrobialMetabolism(targetsFileNameInSDF, preprocess, 
				filter, nr_of_steps, scoreThreshold);		
		this.saveBioTransformationProductsToSdf(biotransformations, metabolitesFileNameInSDF, annotate);
		
	}

	public void  simulateGutMicrobialMetabolismAndSaveToSDF(IAtomContainerSet containers, int nrOfSteps, Double scoreThreshold, String outputFolder) throws Exception {
		simulateGutMicrobialMetabolismAndSaveToSDF(containers, nrOfSteps, scoreThreshold, outputFolder, false);
	}
	
	
	public void  simulateGutMicrobialMetabolismAndSaveToSDF(IAtomContainerSet containers, int nrOfSteps, Double scoreThreshold, String outputFolder, boolean annotate) throws Exception {
		
		if(!containers.isEmpty()){
			for(IAtomContainer molecule : containers.atomContainers()){

				String identifier = molecule.getProperty(CDKConstants.TITLE);
				if(identifier == null){
					identifier = molecule.getProperty("Name");
					if(identifier == null){
						identifier = molecule.getProperty("InChiKey");
						if(identifier == null){
							identifier = this.inchiGenFactory.getInChIGenerator(molecule).getInchiKey();
						}
					}
				}
				identifier = identifier.replace(":", "-").replace("/", "_");
				System.out.println(identifier);
				ArrayList<Biotransformation> biotransformations = this.simulateGutMicrobialMetabolism(molecule, true, true, nrOfSteps, scoreThreshold);
				System.out.println(biotransformations.size() + " biotransformations.");
				this.saveBioTransformationProductsToSdf(biotransformations, outputFolder + "/" + identifier + "_EC_based_metabolites.sdf", annotate);
			
			}
		}		
	}
	
	
	
	public boolean isDeconjugationCandidate(IAtomContainer molecule) throws SMARTSException, CDKException, IOException{
		boolean dec = false;
		
		for(MetabolicReaction m : this.reactionsByGroups.get("deconjugationReactions")){
			if(ChemStructureExplorer.compoundMatchesReactionConstraints(m, molecule)){
				dec = true;
				break;
			}
		}
		
		return dec;		
	}

	public boolean isConjugationCandidate(IAtomContainer molecule) throws SMARTSException, CDKException, IOException{
		boolean dec = false;
		
		for(MetabolicReaction m : this.reactionsByGroups.get("gutMicroPhaseIIReactions")){
			if(ChemStructureExplorer.compoundMatchesReactionConstraints(m, molecule)){
				dec = true;
				break;
			}
		}
		
		return dec;		
	}
	
	
	public void printStatistics(){
		int count = 0;
		for(Enzyme e: this.enzymesList){
			count = count + e.getReactionsNames().size();
		}
		System.out.println("Humber of enzymes: " + this.enzymesList.size());
		System.out.println("Humber of biotransformation rules: " + (this.reactionsByGroups.get("gutMicroReductiveReactions").size() 
				+ this.reactionsByGroups.get("gutMicroPhaseIIReactions").size() 
				+ this.reactionsByGroups.get("deconjugationReactions").size()));		
		System.out.println("Humber of enzyme-biotransformation rules associations: " + count);	
	}
}
