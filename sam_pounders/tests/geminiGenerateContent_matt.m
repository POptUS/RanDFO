function response = geminiGenerateContent_matt(prompt, my_json)

    if strlength(my_json) == 0
        model = "gemini-2.5-flash"; %"gemini-pro";
        query = struct("contents",[]);
        query.contents = {struct("parts",[])};
        query.contents{1}.parts{1} = {struct("text",prompt)};
    else
        model = "gemini-2.5-flash";
        % fid = fopen(img);
        % im = fread(fid,'*uint8');
        % fclose(fid);
        % b64 = matlab.net.base64encode(im);
        % [~,~,ext] = fileparts(img);
        % MIMEType = "image/" + erase(ext,".");
        MIMEType = "application/json";
        query = struct("contents",[]);
        query.contents = {struct("parts",[])};
        query.contents{1}.parts = {struct("text",prompt),struct("inline_data",[])};
        query.contents{1}.parts{2}.inline_data = struct("mime_type",MIMEType,"data",[]);
        query.contents{1}.parts{2}.inline_data.data = my_json;
        % if isfile("imdata.png")
        %     delete("imdata.png")
        % end
    end
    endpoint = "https://generativelanguage.googleapis.com/v1beta/";
    method = "generateContent";
    
    import matlab.net.*
    import matlab.net.http.*  
    apikey = 'AIzaSyC72LqqXVv1auYzhFh5DRemvO7nzqD6v5Q'; %getenv("GEMINI_API_KEY"); 
    headers = HeaderField('Content-Type', 'application/json');
    request = RequestMessage('post', headers, query);
    response = send(request, URI(endpoint + ...
        "models/" + model + ...
        ":" + method + ...
        "?key=" + apikey));
end