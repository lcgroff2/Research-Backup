/**
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */

package biotransformer.biosystems;

import java.io.IOException;

import org.codehaus.jackson.map.ObjectMapper;
import org.json.simple.parser.ParseException;

import biotransformer.biosystems.BioSystem.BioSystemName;

public class Human extends BioSystem {

	public Human(ObjectMapper mapper) throws IOException, ParseException {
		super(BioSystemName.HUMAN, mapper);
	}

}
