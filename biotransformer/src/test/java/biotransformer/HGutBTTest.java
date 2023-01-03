package biotransformer;

import java.io.IOException;
import java.util.ArrayList;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import ambit2.smarts.SMIRKSReaction;
import biotransformer.biomolecule.Enzyme;
import biotransformer.btransformers.HGutBTransformer;
import biotransformer.transformation.Biotransformation;
import biotransformer.transformation.MReactionSets;
import biotransformer.transformation.MetabolicReaction;
import biotransformer.transformation.MRPatterns.ReactionName;
import biotransformer.utils.ChemStructureExplorer;

public class HGutBTTest extends HGutBTransformer{

	public HGutBTTest() throws IOException, CDKException{
		// TODO Auto-generated constructor stub
	}

	
	public static void main(String[] args) throws Exception{
		HGutBTTest em = new HGutBTTest();
		System.out.println(em.getBioSystemName());
		IAtomContainer mol2 = em.getSmiParser().parseSmiles("O[C@@H]1CC2=C(O)C=C(O)C=C2O[C@@H]1C1=CC=C(O)C(O)=C1");
		ArrayList<Biotransformation> biotransformations = em.simulateGutMicrobialMetabolism(mol2, true, true, 2, 0.0);
//		IAtomContainerSet acMetabolites = em.extractAtomContainer(biotransformations);
		
		for(Biotransformation b : biotransformations){
			b.display();
			System.out.println();
		}

	}	
	


}
