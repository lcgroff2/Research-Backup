/**
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */

package biotransformer.utils;

import java.util.ArrayList;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import ambit2.smarts.SMIRKSManager;
import biotransformer.transformation.MReactionSets;
import biotransformer.transformation.MetabolicReaction;

public class ChemStructureManipulator {
	protected static SMIRKSManager smrkMan = new SMIRKSManager(SilentChemObjectBuilder.getInstance());

	
	public ChemStructureManipulator() {
		// TODO Auto-generated constructor stub
		smrkMan.setFlagApplyStereoTransformation(false);
		smrkMan.setFlagCheckResultStereo(true);
		smrkMan.setFlagFilterEquivalentMappings(true);
		smrkMan.setFlagProcessResultStructures(true);
		smrkMan.setFlagAddImplicitHAtomsOnResultProcess(true);
	}
	
	/**
	 * This function applies some preprocessing operations, such as setting the
	 * flag of atoms from aromatic rings to "ISAROMATIC", and kelulizing
	 * molecules.
	 * 
	 * @param molecule
	 *            : A molecule of interest
	 * @return : A processed molecule (AtomContainer)
	 * @throws CDKException
	 */
	public static IAtomContainer preprocessContainer(IAtomContainer molecule)
			throws CDKException {
		
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
//		AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
		
//	    Aromaticity aromaticity = new Aromaticity( ElectronDonation.cdk(), Cycles.cdkAromaticSet());
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.daylight(), Cycles.or(Cycles.all(), Cycles.all(6)));
		
		for (IBond bond : molecule.bonds()) {
			if (bond.isAromatic() && bond.getOrder() == IBond.Order.UNSET) {
				bond.setFlag(CDKConstants.ISAROMATIC, true);
				bond.getAtom(0).setFlag(CDKConstants.ISAROMATIC, true);
				bond.getAtom(1).setFlag(CDKConstants.ISAROMATIC, true);

			} 
		}
//		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		aromaticity.apply(molecule);
//		Kekulization.kekulize(molecule);
//		AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
		StructureDiagramGenerator sdg = new StructureDiagramGenerator();
		sdg.setMolecule(molecule);
		sdg.generateCoordinates();		
		IAtomContainer layedOutMol = sdg.getMolecule();

//		StringWriter w2 = new StringWriter();
//		MDLWriter mw2 = new MDLWriter(w2);
//		mw2.write(layedOutMol);	
//		System.out.println("After preprocessing\n" + w2.toString() + "\n\n");

		return layedOutMol;
	}
	
	public static IAtomContainer standardizeMoleculeWithCopy(IAtomContainer molecule) throws Exception{
		return  standardizeMoleculeWithCopy(molecule, true);
	}


	public static IAtomContainer standardizeMoleculeWithCopy(IAtomContainer molecule, boolean preprocess) throws Exception{
	//	System.out.println("Starting the standardization process...");
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance());
		IAtomContainer molClone =  molecule.clone();
		if(preprocess){
			molClone = ChemStructureManipulator.preprocessContainer(molClone);
	//		System.out.println("INCHIKEY BEFORE ADDING EXPLICIT HYDROGEN: " +  this.inchiGenFactory.getInChIGenerator(molClone).getInchiKey());
	//		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molClone);
			AtomContainerManipulator.convertImplicitToExplicitHydrogens(molClone);
	//		System.out.println("INCHIKEY AFTER ADDING EXPLICIT HYDROGEN: " +  this.inchiGenFactory.getInChIGenerator(molClone).getInchiKey());
		} else {
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molClone);
			AtomContainerManipulator.convertImplicitToExplicitHydrogens(molClone);
		}
		
		ArrayList<MetabolicReaction> matchedReactions = new ArrayList<MetabolicReaction>();
		for(MetabolicReaction mr : MReactionSets.standardizationReactions){
			if(ChemStructureExplorer.compoundMatchesReactionConstraints(mr, molClone)){
				matchedReactions.add(mr);
			}
		}
		
	//	molClone = AtomContainerManipulator.suppressHydrogens(molClone);
	//	System.out.println("The processed molecule before standardization is: " + this.smiGen.isomeric().create(molClone));
		
		for(MetabolicReaction m : matchedReactions){
	//		System.out.println(m.name);
			while(ChemStructureExplorer.compoundMatchesReactionConstraints(m, molClone)){
//				System.out.println("Molecule " + ChemStructureExplorer.smiGen.create(molClone) + " still matches: " + m.name);
	//			molClone = this.smrkMan.applyTransformationWithSingleCopyForEachPos(molClone, null, m.getSmirksReaction()).getAtomContainer(0);
				smrkMan.applyTransformation(molClone, m.getSmirksReaction());
				AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molClone);
				
				// I added the following line because NullPointerExceptions were returned complaining about some atoms with unset implicit hydrogens
				// This occured when the comppound was transfromed.
				
				adder.addImplicitHydrogens(molClone);
//				System.out.println("The processed molecule after " + m.getReactionName() + " is : " + ChemStructureExplorer.smiGen.create(molClone));
			
			}

		}

		AtomContainerManipulator.suppressHydrogens(molClone);
		return molClone;
	}
	
	public static IAtomContainer removeHydrogenAtoms(IAtomContainer atc){
		AtomContainerManipulator.removeHydrogens(atc);
		return atc;
	}
	
	
//	public static boolean isUnwantedCompound(IAtomContainer mol){
//		boolean unwantedCompound = false;
//		
//		
//		return unwantedCompound;
//		
//	}
	
	

}
