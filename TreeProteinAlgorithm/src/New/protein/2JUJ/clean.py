with open("2JUJ.pdb") as f_in:
		with open("2JUJ.clean.pdb","w") as f_out:
			for line in f_in:
				if line[:4]=="ATOM":
					f_out.write(line)
