package biotransformer;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.Time;
import java.util.ArrayList;

import org.json.JSONObject;
import org.json.simple.parser.ParseException;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import com.mashape.unirest.http.HttpResponse;
import com.mashape.unirest.http.JsonNode;
import com.mashape.unirest.http.Unirest;

import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.btransformers.Biotransformer;
import biotransformer.transformation.Biotransformation;
import biotransformer.transformation.MetabolicReaction;
import biotransformer.transformation.MRPatterns.ReactionName;
import biotransformer.utils.ChemStructureExplorer;
import biotransformer.utils.ChemStructureManipulator;
import biotransformer.utils.ChemdbRest;
import biotransformer.utils.FileUtilities;
import biotransformer.utils.HumanSuperBioTransformer;


public class HumanSuperBTransformerTest extends HumanSuperBioTransformer{

	public HumanSuperBTransformerTest() throws IOException, ParseException, CDKException {
		// TODO Auto-generated constructor stub
		super();
	}
	
	public static void main(String[] args) throws Exception{
		HumanSuperBTransformerTest hsbt 	= new HumanSuperBTransformerTest();
		IChemObjectBuilder 	builder 	= SilentChemObjectBuilder.getInstance();
		SmilesGenerator smiGen 			= new SmilesGenerator().isomeric();
		SmilesParser	smiParser		= new SmilesParser(builder);
		

		IAtomContainer molecule = smiParser.parseSmiles("O[C@@H]1CC2=C(O)C=C(O)C=C2O[C@@H]1C1=CC=C(O)C(O)=C1");
		Biotransformer b = new Biotransformer(BioSystemName.HUMAN);
		IAtomContainer mt = ChemStructureManipulator.standardizeMoleculeWithCopy(molecule);
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(mt);
		hsbt.simulateHumanSuperbioMetabolismAndSaveToCSV(molecule, "data/epicatechin-superbio-2.csv", false);
//		hsbt.simulateHumanAndGutMicrobialMetabolismAndSaveToSDF(molecule, "data/epicatechin-superbio-2.sdf", false);
				
	}

}
