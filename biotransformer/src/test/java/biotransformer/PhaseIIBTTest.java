package biotransformer;

import java.io.IOException;
import java.util.ArrayList;

import org.json.simple.parser.ParseException;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;

import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.btransformers.Phase2BTransformer;
import biotransformer.transformation.Biotransformation;
import biotransformer.transformation.MetabolicReaction;
import biotransformer.utils.ChemStructureExplorer;
import biotransformer.utils.ChemStructureManipulator;
import biotransformer.utils.Utilities;
import predicition.P2Filter;

public class PhaseIIBTTest extends Phase2BTransformer{

	public PhaseIIBTTest(BioSystemName bioSName) throws IOException, ParseException, CDKException{
		// TODO Auto-generated constructor stub
		super(bioSName);
	}

		
	public static void main(String[] args) throws Exception{
		PhaseIIBTTest em = new PhaseIIBTTest(BioSystemName.HUMAN);
		P2Filter p2f = new P2Filter();		
		IAtomContainer mol = em.getSmiParser().parseSmiles("O[C@@H]1CC2=C(O)C=C(O)C=C2O[C@@H]1C1=CC=C(O)C(O)=C1");
		IAtomContainer stac = ChemStructureManipulator.standardizeMoleculeWithCopy(mol, true);		
		ArrayList<Biotransformation> biotransformations = em.applyPhase2TransformationsChainAndReturnBiotransformations(mol, true, true, true, 1, 0.0);
		System.out.println("Number of biotransformations: " + biotransformations.size());
		em.saveBioTransformationProductsToSdf(biotransformations, "data/epicatechin_PhaseII_metabolites.sdf", false);	
	}

}
