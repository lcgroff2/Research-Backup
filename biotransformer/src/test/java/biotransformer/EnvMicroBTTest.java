package biotransformer;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import org.json.simple.parser.ParseException;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;

import biotransformer.btransformers.EnvMicroBTransformer;
import biotransformer.transformation.Biotransformation;
import biotransformer.utils.ChemStructureManipulator;
import biotransformer.utils.FileUtilities;

public class EnvMicroBTTest extends EnvMicroBTransformer{

	public EnvMicroBTTest() throws IOException, ParseException, CDKException{
		// TODO Auto-generated constructor stub
		super();
	}


	public static void main(String[] args) throws Exception{
		EnvMicroBTTest em = new EnvMicroBTTest();
		System.out.println(em.getBioSystemName());
		System.out.println(em.getReactionsList().size());
		
		
		IAtomContainer mol = em.getSmiParser().parseSmiles("O[C@@H]1CC2=C(O)C=C(O)C=C2O[C@@H]1C1=CC=C(O)C(O)=C1");
		IAtomContainer stac = ChemStructureManipulator.standardizeMoleculeWithCopy(mol, true);
		IAtomContainerSet filteredCompounds = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		filteredCompounds.addAtomContainer(mol);
		em.simulateEnvMicrobialDegradationAndSaveToSDF(stac, true, true, 2, 1.0, "data", true);

	}
	
}
