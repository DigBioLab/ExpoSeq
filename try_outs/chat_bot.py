import openai



def askMe(global_params):
    openai.api_key = global_params["api_gpt3"]
    if openai.api_key == "":
        api = input("enter your api key from openai. Therefore visit: https://openai.com/")
        global_params["api_gpt3"] = api
        with open("../src/expoSeq/settings/global_vars.txt", "w") as f:
            f.write(str(global_params))
        openai.api_key = api
    else:
        pass
    introduction = openai.Completion.create(model="text-davinci-003",
                             prompt="Welcome the user and tell him or her how you can help out", temperature=0,
                             max_tokens=1000)
    prompt = input(introduction.choices[0].text)
    while True:
        response = openai.Completion.create(model="text-davinci-003", prompt=prompt, temperature=0, max_tokens=1000)
        print(response.choices[0].text)
        conversation= openai.Completion.create(model="text-davinci-003",
                                                prompt="Ask the user if he or she has further questions. Tell the user that he should press enter or type no if he wants to leave the chat",
                                                temperature=0,
                                                max_tokens=1000)
        print(conversation.choices[0].text)
        answer = input()
        if answer.lower() in ["", "No", "no"]:
            break
        else:
            continue
