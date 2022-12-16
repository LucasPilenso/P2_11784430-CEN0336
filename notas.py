total = 0
contador= 0

print("Insira 10 notas entre 0 e 10.") #Informar o usuário das restrições 
                                       #das entradas
while contador < 10:
    try:
        nota = float(input("Insira a nota: "))
        if nota >= 0 and nota <= 10:
            TOTAL += nota              #Gerar o loop dentro do bloco try
            CONTADOR_NOTAS += 1
        else:
            print("Nota inválida. Insira uma nota entre 0 e 10.")
    except (ValueError, TypeError):
        print("Insira um valor válido para a nota.")
                                      #Esclarecer os erros(value e type)
media_disciplina = TOTAL / 10
print("A media da disciplina é: ", media_disciplina)
                                      #Em caso de sucesso, retorna a média
